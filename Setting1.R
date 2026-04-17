# =============================================================================
# Setting 1: Validation of Auto-g-Computation
# =============================================================================
#
# PURPOSE
# -------
# This script validates the auto-g-computation estimator for causal estimands
# (Direct Effect, Indirect Effect, Average Total Effect) in a network setting
# with binary treatment, covariate, and outcome. For each of three graph
# specifications (a), (b), (c), the simulation:
#
#   1. Computes the analytical (ground-truth) DE, IE, ATE via enumeration
#      (feasible because groups have only 4 units).
#   2. Runs n_sim Monte Carlo replicates, where each replicate:
#        a. Draws treatment allocations from the graph's OWN observational
#           Gibbs sampler.
#        b. Runs the interventional Gibbs sampler under the SAME graph.
#        c. Computes DE, IE, ATE estimates.
#   3. Reports empirical mean, bias, empirical SE, RMSE, and 95% coverage.
#
# GRAPHS
# ------
# Each group of 4 units uses one of two local adjacency patterns:
#   - chain : 1 - 2 - 3 - 4
#   - pairs : (1 - 2), (3 - 4)
#
#   Graph (a) : interference = chain, dependence = chain
#   Graph (b) : interference = pairs, dependence = pairs
#   Graph (c) : interference = pairs, dependence = chain
#
# COVERAGE CONVENTION
# -------------------
#   coverage = mean( | est - truth | <= 1.96 * sd(est) )
# where sd(est) is a single across-replicate empirical standard deviation.
#
# NOTE
# ----
# This script uses the same generic parameterized pipeline as the Setting 2
# code: one set of energy / Gibbs functions, parameterized by the local
# interference and dependence adjacency matrices (G_int_local, G_dep_local).
#
# AUTHOR:  <your name>
# DATE  :  <date>
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Dependencies
# -----------------------------------------------------------------------------
library(parallel)
library(doParallel)
library(foreach)


# -----------------------------------------------------------------------------
# 1. Network constructors
# -----------------------------------------------------------------------------
# The population of N units is partitioned into N/4 disjoint groups of size 4.
# The full N x N adjacency matrices are block-diagonal, with each 4 x 4 block
# given by a "local" adjacency pattern (chain or pairs).
# -----------------------------------------------------------------------------

#' Build block-diagonal interference and dependence adjacency matrices.
#'
#' @param N            Total number of units (must be divisible by 4).
#' @param G_int_local  4x4 local interference adjacency.
#' @param G_dep_local  4x4 local dependence adjacency.
#' @return A list with the full N x N matrices plus bookkeeping fields.
build_network <- function(N, G_int_local, G_dep_local) {
  stopifnot(N %% 4 == 0)
  n_groups   <- N / 4
  group_size <- 4
  G_int <- matrix(0, N, N)
  G_dep <- matrix(0, N, N)
  for (g in 1:n_groups) {
    idx <- ((g - 1) * group_size + 1):(g * group_size)
    G_int[idx, idx] <- G_int_local
    G_dep[idx, idx] <- G_dep_local
  }
  list(G_interference = G_int,
       G_dependence   = G_dep,
       N              = N,
       n_groups       = n_groups,
       group_size     = group_size)
}

#' Chain adjacency 1 - 2 - 3 - 4 (used as dependence in graph c; both in a).
chain_local <- function() {
  G <- matrix(0, 4, 4)
  for (k in 1:3) G[k, k + 1] <- G[k + 1, k] <- 1
  G
}

#' Pairs adjacency (1-2), (3-4) (used as interference in graph c; both in b).
pairs_local <- function() {
  G <- matrix(0, 4, 4)
  G[1, 2] <- G[2, 1] <- 1
  G[3, 4] <- G[4, 3] <- 1
  G
}

# Convenience constructors for the three graph specifications.
create_network_a <- function(N = 8) build_network(N, chain_local(), chain_local())
create_network_b <- function(N = 8) build_network(N, pairs_local(), pairs_local())
create_network_c <- function(N = 8) build_network(N, pairs_local(), chain_local())


# -----------------------------------------------------------------------------
# 2. Utilities
# -----------------------------------------------------------------------------

#' Enumerate all 2^n binary vectors as a 2^n x n matrix.
generate_binary_configs <- function(n) {
  configs <- matrix(0, 2^n, n)
  for (j in 1:n)
    configs[, j] <- rep(c(0, 1), each = 2^(n - j), times = 2^(j - 1))
  configs
}

#' Degree of node i in a 4x4 adjacency matrix, clamped to >= 1 so that
#' reciprocal-degree weights never divide by zero.
deg_local <- function(G_local, i) {
  d <- sum(G_local[i, ])
  if (d == 0) 1 else d
}


# -----------------------------------------------------------------------------
# 3. Group-level energy functions
# -----------------------------------------------------------------------------
# These define the joint distributions of (L, A, Y) within a single group of 4.
#
# Parameterization:
#   - Directed contributions (neighbor covariates in L, A): mean() over
#     interference neighbors.
#   - Undirected contributions (Y-Y, A-A, L-L): symmetric weight
#     w_avg = 0.5 * (1/deg_i + 1/deg_j).
# -----------------------------------------------------------------------------

#' Group energy for outcome Y given (A, L) within a group.
group_energy_Y <- function(Y_vec, A_vec, L_vec, params,
                           G_int_local, G_dep_local) {
  energy <- 0
  # Node-level linear predictor terms
  for (i in 1:4) {
    nb_int <- which(G_int_local[i, ] == 1)
    lp <- params$beta0 + params$beta1 * A_vec[i] + params$beta2 * L_vec[i]
    if (length(nb_int) > 0) {
      lp <- lp +
        params$beta3 * mean(A_vec[nb_int]) +
        params$beta4 * mean(L_vec[nb_int])
    }
    energy <- energy + Y_vec[i] * lp
  }
  # Pairwise Y-Y dependence terms
  for (i in 1:3) {
    for (j in (i + 1):4) {
      if (G_dep_local[i, j] == 1) {
        w_avg <- 0.5 * (1 / deg_local(G_dep_local, i) +
                        1 / deg_local(G_dep_local, j))
        energy <- energy + Y_vec[i] * Y_vec[j] * w_avg * params$theta
      }
    }
  }
  energy
}

#' Group energy for treatment A given L within a group.
group_energy_A <- function(A_vec, L_vec, params_A,
                           G_int_local, G_dep_local) {
  energy <- 0
  for (i in 1:4) {
    nb_int <- which(G_int_local[i, ] == 1)
    energy <- energy + A_vec[i] * (params_A$gamma0 + params_A$gamma1 * L_vec[i])
    if (length(nb_int) > 0) {
      energy <- energy + A_vec[i] * params_A$gamma2 * mean(L_vec[nb_int])
    }
  }
  for (i in 1:3) {
    for (j in (i + 1):4) {
      if (G_dep_local[i, j] == 1) {
        w_avg <- 0.5 * (1 / deg_local(G_dep_local, i) +
                        1 / deg_local(G_dep_local, j))
        energy <- energy + A_vec[i] * A_vec[j] * w_avg * params_A$psi
      }
    }
  }
  energy
}

#' Group energy for covariate L within a group.
group_energy_L <- function(L_vec, params_L, G_dep_local) {
  energy <- sum(L_vec) * params_L$alpha
  for (i in 1:3) {
    for (j in (i + 1):4) {
      if (G_dep_local[i, j] == 1) {
        w_avg <- 0.5 * (1 / deg_local(G_dep_local, i) +
                        1 / deg_local(G_dep_local, j))
        energy <- energy + L_vec[i] * L_vec[j] * w_avg * params_L$omega
      }
    }
  }
  energy
}


# -----------------------------------------------------------------------------
# 4. Full conditional probabilities (for Gibbs sampling)
# -----------------------------------------------------------------------------
# Each full conditional P(X_i = 1 | rest) is obtained by computing the group
# energy under X_i = 1 vs X_i = 0, then returning a numerically stable softmax.
# -----------------------------------------------------------------------------

#' P(Y_i = 1 | Y_{-i}, A, L).
prob_Y_given_rest <- function(i, Y, A, L, params, network,
                              G_int_local, G_dep_local) {
  gs    <- network$group_size
  g     <- ceiling(i / gs)
  i_loc <- i - (g - 1) * gs
  idx   <- ((g - 1) * gs + 1):(g * gs)
  Yg <- Y[idx]; Ag <- A[idx]; Lg <- L[idx]
  Y1 <- Yg; Y1[i_loc] <- 1
  Y0 <- Yg; Y0[i_loc] <- 0
  e1 <- group_energy_Y(Y1, Ag, Lg, params, G_int_local, G_dep_local)
  e0 <- group_energy_Y(Y0, Ag, Lg, params, G_int_local, G_dep_local)
  me <- max(e0, e1)                          # numerical stabilizer
  exp(e1 - me) / (exp(e0 - me) + exp(e1 - me))
}

#' P(A_i = 1 | A_{-i}, L).
prob_A_given_rest <- function(i, A, L, params_A, network,
                              G_int_local, G_dep_local) {
  gs    <- network$group_size
  g     <- ceiling(i / gs)
  i_loc <- i - (g - 1) * gs
  idx   <- ((g - 1) * gs + 1):(g * gs)
  Ag <- A[idx]; Lg <- L[idx]
  A1 <- Ag; A1[i_loc] <- 1
  A0 <- Ag; A0[i_loc] <- 0
  e1 <- group_energy_A(A1, Lg, params_A, G_int_local, G_dep_local)
  e0 <- group_energy_A(A0, Lg, params_A, G_int_local, G_dep_local)
  me <- max(e0, e1)
  exp(e1 - me) / (exp(e0 - me) + exp(e1 - me))
}

#' P(L_i = 1 | L_{-i}).
prob_L_given_rest <- function(i, L, params_L, network, G_dep_local) {
  gs    <- network$group_size
  g     <- ceiling(i / gs)
  i_loc <- i - (g - 1) * gs
  idx   <- ((g - 1) * gs + 1):(g * gs)
  Lg <- L[idx]
  L1 <- Lg; L1[i_loc] <- 1
  L0 <- Lg; L0[i_loc] <- 0
  e1 <- group_energy_L(L1, params_L, G_dep_local)
  e0 <- group_energy_L(L0, params_L, G_dep_local)
  me <- max(e0, e1)
  exp(e1 - me) / (exp(e0 - me) + exp(e1 - me))
}


# -----------------------------------------------------------------------------
# 5. Analytical ground truth
# -----------------------------------------------------------------------------
# Because each group has only 4 units, we can enumerate all 2^4 = 16 binary
# configurations and compute the true DE, IE, and ATE in closed form.
# -----------------------------------------------------------------------------

#' Exact E[Y_i | do(A = a_vec)] by marginalizing L and Y analytically.
psi_analytical_group <- function(i_loc, a_vec, params, params_L,
                                 G_int_local, G_dep_local) {
  L_cfg <- generate_binary_configs(4)
  Y_cfg <- generate_binary_configs(4)
  # Marginal distribution of L
  fL <- exp(apply(L_cfg, 1, function(l)
    group_energy_L(l, params_L, G_dep_local)))
  fL <- fL / sum(fL)
  # E[Y_i | A = a_vec, L = l] for each l
  EYi <- sapply(1:nrow(L_cfg), function(li) {
    fY <- exp(apply(Y_cfg, 1, function(y)
      group_energy_Y(y, a_vec, L_cfg[li, ], params, G_int_local, G_dep_local)))
    sum(Y_cfg[, i_loc] * fY / sum(fY))
  })
  sum(EYi * fL)
}

#' Exact P(A = A_vec | L = L_vec).
prob_A_given_L_analytical <- function(A_vec, L_vec, params_A,
                                      G_int_local, G_dep_local) {
  A_cfg <- generate_binary_configs(4)
  fA <- exp(apply(A_cfg, 1, function(a)
    group_energy_A(a, L_vec, params_A, G_int_local, G_dep_local)))
  a_idx <- which(apply(A_cfg, 1, function(x) all(x == A_vec)))
  fA[a_idx] / sum(fA)
}

#' Marginal allocation probability pi(a) = P(A = a_vec) under the observational
#' data-generating process.
marginal_pi_A <- function(A_vec, params_A, params_L,
                          G_int_local, G_dep_local) {
  L_cfg <- generate_binary_configs(4)
  fL <- exp(apply(L_cfg, 1, function(l)
    group_energy_L(l, params_L, G_dep_local)))
  fL <- fL / sum(fL)
  sum(sapply(1:nrow(L_cfg), function(li)
    prob_A_given_L_analytical(A_vec, L_cfg[li, ], params_A,
                              G_int_local, G_dep_local) * fL[li]))
}

#' Compute the analytical (ground-truth) DE, IE, and ATE for a given graph.
calculate_estimands_analytical <- function(network, params, params_A, params_L,
                                           G_int_local, G_dep_local,
                                           verbose = TRUE) {
  gs    <- network$group_size
  A_cfg <- generate_binary_configs(gs)

  if (verbose) cat("  Computing allocation probabilities pi(a)...\n")
  pi_A <- sapply(1:nrow(A_cfg), function(k)
    marginal_pi_A(A_cfg[k, ], params_A, params_L, G_int_local, G_dep_local))
  if (verbose) cat("  Allocation probabilities sum to:", round(sum(pi_A), 6), "\n")

  if (verbose) cat("  Computing analytical psi for each position...\n")
  psi_1 <- psi_0 <- psi_0_0 <- numeric(gs)
  for (i in 1:gs) {
    if (verbose) cat("    Position", i, "\n")
    # psi(1, i): average over allocations with A_i = 1
    cfg1 <- which(A_cfg[, i] == 1); pi1 <- pi_A[cfg1] / sum(pi_A[cfg1])
    psi_1[i] <- sum(pi1 * sapply(cfg1, function(k)
      psi_analytical_group(i, A_cfg[k, ], params, params_L,
                           G_int_local, G_dep_local)))
    # psi(0, i): average over allocations with A_i = 0
    cfg0 <- which(A_cfg[, i] == 0); pi0 <- pi_A[cfg0] / sum(pi_A[cfg0])
    psi_0[i] <- sum(pi0 * sapply(cfg0, function(k)
      psi_analytical_group(i, A_cfg[k, ], params, params_L,
                           G_int_local, G_dep_local)))
    # psi(0, 0): everyone untreated
    psi_0_0[i] <- psi_analytical_group(i, rep(0, gs), params, params_L,
                                       G_int_local, G_dep_local)
  }

  DE  <- mean(psi_1 - psi_0)      # Direct effect
  IE  <- mean(psi_0 - psi_0_0)    # Indirect (spillover) effect
  ATE <- mean(psi_1 - psi_0_0)    # Average total effect

  list(DE = DE, IE = IE, ATE = ATE)
}


# -----------------------------------------------------------------------------
# 6. Gibbs samplers
# -----------------------------------------------------------------------------
# Used for (i) generating observational data from which to draw allocations,
# and (ii) simulating Y | do(A = a_vec) for a chosen allocation.
# -----------------------------------------------------------------------------

#' Observational Gibbs sampler for (L, A, Y).
gibbs_sampler_observational <- function(network, params, params_A, params_L,
                                        G_int_local, G_dep_local,
                                        n_iter = 5000, burn_in = 1000) {
  N <- network$N
  L <- rbinom(N, 1, 0.5); A <- rbinom(N, 1, 0.5); Y <- rbinom(N, 1, 0.5)
  n_s <- n_iter - burn_in
  Ls <- matrix(0, n_s, N); As <- matrix(0, n_s, N); Ys <- matrix(0, n_s, N)
  si <- 1
  for (m in 1:n_iter) {
    for (i in sample(1:N)) {
      L[i] <- rbinom(1, 1, prob_L_given_rest(i, L, params_L, network, G_dep_local))
      A[i] <- rbinom(1, 1, prob_A_given_rest(i, A, L, params_A, network,
                                             G_int_local, G_dep_local))
      Y[i] <- rbinom(1, 1, prob_Y_given_rest(i, Y, A, L, params, network,
                                             G_int_local, G_dep_local))
    }
    if (m > burn_in) { Ls[si, ] <- L; As[si, ] <- A; Ys[si, ] <- Y; si <- si + 1 }
  }
  list(L_samples = Ls, A_samples = As, Y_samples = Ys)
}

#' Interventional Gibbs sampler for (L, Y) with A held fixed at a_vec.
gibbs_sampler_interventional <- function(a_vec, network, params, params_L,
                                         G_int_local, G_dep_local,
                                         n_iter = 5000, burn_in = 1000) {
  N <- network$N
  L <- rbinom(N, 1, 0.5); Y <- rbinom(N, 1, 0.5)
  n_s <- n_iter - burn_in
  Ls <- matrix(0, n_s, N); Ys <- matrix(0, n_s, N)
  si <- 1
  for (m in 1:n_iter) {
    for (i in sample(1:N)) {
      L[i] <- rbinom(1, 1, prob_L_given_rest(i, L, params_L, network, G_dep_local))
      Y[i] <- rbinom(1, 1, prob_Y_given_rest(i, Y, a_vec, L, params, network,
                                             G_int_local, G_dep_local))
    }
    if (m > burn_in) { Ls[si, ] <- L; Ys[si, ] <- Y; si <- si + 1 }
  }
  list(L_samples = Ls, Y_samples = Ys)
}


# -----------------------------------------------------------------------------
# 7. Single replicate
# -----------------------------------------------------------------------------
# Workflow for one Monte-Carlo replicate on a single graph:
#   (1) Draw n_alloc allocations a_vec from THAT graph's observational sampler.
#   (2) For each allocation, run an interventional sampler under the SAME graph.
#   (3) Aggregate into DE / IE / ATE estimates.
# -----------------------------------------------------------------------------
run_single_replicate <- function(sim_id,
                                 network, N,
                                 G_int_local, G_dep_local,
                                 params, params_A, params_L,
                                 n_iter, burn_in, n_alloc,
                                 seed = NULL) {
  if (!is.null(seed)) set.seed(seed + sim_id)

  # Step 1: draw allocations from the graph's own observational sampler
  obs <- gibbs_sampler_observational(network, params, params_A, params_L,
                                     G_int_local, G_dep_local,
                                     n_iter = n_alloc * 50, burn_in = 500)
  aidx <- seq(1, nrow(obs$A_samples), length.out = n_alloc)

  # Step 2: interventional estimation for each allocation
  psi_1 <- psi_0 <- cnt1 <- cnt0 <- numeric(N)
  for (k in 1:n_alloc) {
    a_vec <- obs$A_samples[round(aidx[k]), ]
    Ymeans <- colMeans(
      gibbs_sampler_interventional(a_vec, network, params, params_L,
                                   G_int_local, G_dep_local,
                                   n_iter, burn_in)$Y_samples
    )
    for (i in 1:N) {
      if (a_vec[i] == 1) {
        psi_1[i] <- psi_1[i] + Ymeans[i]; cnt1[i] <- cnt1[i] + 1
      } else {
        psi_0[i] <- psi_0[i] + Ymeans[i]; cnt0[i] <- cnt0[i] + 1
      }
    }
  }
  psi_1 <- psi_1 / pmax(cnt1, 1)
  psi_0 <- psi_0 / pmax(cnt0, 1)

  # All-control counterfactual
  psi_0_0 <- colMeans(
    gibbs_sampler_interventional(rep(0, N), network, params, params_L,
                                 G_int_local, G_dep_local,
                                 n_iter, burn_in)$Y_samples
  )

  c(DE  = mean(psi_1 - psi_0),
    IE  = mean(psi_0 - psi_0_0),
    ATE = mean(psi_1 - psi_0_0))
}


# -----------------------------------------------------------------------------
# 8. Parallel simulation for a single graph
# -----------------------------------------------------------------------------
# For a given graph specification (G_int_local, G_dep_local):
#   1. Compute the analytical truth.
#   2. Run n_sim replicates in parallel via foreach + doParallel.
#   3. Return truth, raw per-replicate estimates, and elapsed time.
# -----------------------------------------------------------------------------
run_graph_simulation <- function(graph_label,
                                 network, G_int_local, G_dep_local,
                                 n_sim, params, params_A, params_L,
                                 n_iter, burn_in, n_alloc,
                                 n_cores, seed, verbose = TRUE) {
  if (verbose) {
    cat("===========================================\n")
    cat("Graph", graph_label, "| N =", network$N,
        "| n_sim =", n_sim, "| n_cores =", n_cores, "\n")
    cat("===========================================\n\n")
    cat("Calculating analytical truth...\n")
  }
  truth <- calculate_estimands_analytical(network, params, params_A, params_L,
                                          G_int_local, G_dep_local,
                                          verbose = verbose)
  if (verbose) {
    cat("\nTrue DE:",  round(truth$DE, 4),
        "| True IE:", round(truth$IE, 4),
        "| True ATE:", round(truth$ATE, 4), "\n\n")
  }

  # --- Parallel backend -------------------------------------------------------
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  clusterExport(cl, c(
    "build_network", "chain_local", "pairs_local",
    "create_network_a", "create_network_b", "create_network_c",
    "generate_binary_configs", "deg_local",
    "group_energy_Y", "group_energy_A", "group_energy_L",
    "prob_Y_given_rest", "prob_A_given_rest", "prob_L_given_rest",
    "gibbs_sampler_observational", "gibbs_sampler_interventional",
    "run_single_replicate"
  ), envir = environment())

  # --- Monte Carlo loop -------------------------------------------------------
  if (verbose) cat("Running", n_sim, "replicates in parallel...\n")
  t0 <- Sys.time()
  res <- foreach(sim = 1:n_sim, .combine = rbind) %dopar% {
    run_single_replicate(
      sim_id = sim,
      network = network, N = network$N,
      G_int_local = G_int_local, G_dep_local = G_dep_local,
      params = params, params_A = params_A, params_L = params_L,
      n_iter = n_iter, burn_in = burn_in, n_alloc = n_alloc,
      seed = seed
    )
  }
  stopCluster(cl)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  if (verbose) cat("Completed in", round(elapsed, 2), "minutes\n\n")

  colnames(res) <- c("DE", "IE", "ATE")
  list(truth = truth, raw = res, elapsed_mins = elapsed)
}


# -----------------------------------------------------------------------------
# 9. Main driver: run all three graphs and summarize
# -----------------------------------------------------------------------------
# Orchestrates the full Setting 1 simulation by calling run_graph_simulation()
# once per graph (a/b/c), then summarizing bias, empirical SE, RMSE, and 95%
# coverage for each (graph, estimand) combination.
# -----------------------------------------------------------------------------
run_simulation_setting1 <- function(
    n_sim    = 500,
    N        = 8,
    params   = list(beta0 = -1, beta1 = 0.5, beta2 = 0.3,
                    beta3 = 0.4, beta4 = 0.2, theta = 0.3),
    params_A = list(gamma0 = 0, gamma1 = 0.3, gamma2 = 0.2, psi = 0.4),
    params_L = list(alpha = -0.5, omega = 0.4),
    n_iter   = 5000,
    burn_in  = 1000,
    n_alloc  = 50,
    n_cores  = NULL,
    verbose  = TRUE,
    seed     = 42) {

  if (is.null(n_cores)) n_cores <- max(1, detectCores() - 1)

  # --- Graph definitions ------------------------------------------------------
  G_int_a <- chain_local(); G_dep_a <- chain_local()   # graph (a)
  G_int_b <- pairs_local(); G_dep_b <- pairs_local()   # graph (b)
  G_int_c <- pairs_local(); G_dep_c <- chain_local()   # graph (c)

  net_a <- create_network_a(N)
  net_b <- create_network_b(N)
  net_c <- create_network_c(N)

  if (verbose) {
    cat("###############################################\n")
    cat("SETTING 1: Auto-g-computation validation\n")
    cat("###############################################\n\n")
  }

  # --- Run each graph ---------------------------------------------------------
  res_a <- run_graph_simulation("(a)", net_a, G_int_a, G_dep_a,
                                n_sim, params, params_A, params_L,
                                n_iter, burn_in, n_alloc,
                                n_cores, seed, verbose)
  res_b <- run_graph_simulation("(b)", net_b, G_int_b, G_dep_b,
                                n_sim, params, params_A, params_L,
                                n_iter, burn_in, n_alloc,
                                n_cores, seed, verbose)
  res_c <- run_graph_simulation("(c)", net_c, G_int_c, G_dep_c,
                                n_sim, params, params_A, params_L,
                                n_iter, burn_in, n_alloc,
                                n_cores, seed, verbose)

  # --- Summarization ----------------------------------------------------------
  # 95% coverage based on a single across-replicate SD.
  coverage <- function(est, truth_val)
    mean(abs(est - truth_val) <= 1.96 * sd(est))

  summarize_graph <- function(label, result) {
    truth <- result$truth; raw <- result$raw
    do.call(rbind, lapply(c("DE", "IE", "ATE"), function(nm) {
      est <- raw[, nm]; tv <- truth[[nm]]
      data.frame(
        graph    = label,
        estimand = nm,
        truth    = round(tv, 4),
        mean     = round(mean(est), 4),
        bias     = round(mean(est) - tv, 4),
        emp_se   = round(sd(est), 4),
        rmse     = round(sqrt(mean((est - tv)^2)), 4),
        coverage = round(100 * coverage(est, tv), 1),
        stringsAsFactors = FALSE
      )
    }))
  }

  summary_table <- rbind(
    summarize_graph("Graph (a)", res_a),
    summarize_graph("Graph (b)", res_b),
    summarize_graph("Graph (c)", res_c)
  )

  if (verbose) {
    cat("###############################################\n")
    cat("RESULTS - SETTING 1\n")
    cat("###############################################\n")
    print(summary_table, row.names = FALSE)
    cat("\nElapsed per graph (mins):",
        round(res_a$elapsed_mins, 2), "(a) |",
        round(res_b$elapsed_mins, 2), "(b) |",
        round(res_c$elapsed_mins, 2), "(c)\n")
    cat("Total:",
        round(res_a$elapsed_mins + res_b$elapsed_mins + res_c$elapsed_mins, 2),
        "minutes on", n_cores, "cores\n")
  }

  list(graph_a = res_a,
       graph_b = res_b,
       graph_c = res_c,
       summary = summary_table,
       n_cores = n_cores)
}


# =============================================================================
# 10. Entry point
# =============================================================================
# Edit parameter lists below to reproduce / modify the reported run.
# =============================================================================

if (sys.nframe() == 0) {

  params   <- list(beta0 = -1, beta1 = 0.5, beta2 = 0.3,
                   beta3 = 0.4, beta4 = 0.2, theta = 0.3)
  params_A <- list(gamma0 = 0, gamma1 = 0.3, gamma2 = 0.2, psi = 0.4)
  params_L <- list(alpha = -0.5, omega = 0.4)

  cat("Detected", detectCores(), "CPU cores\n\n")

  results_s1 <- run_simulation_setting1(
    n_sim    = 500,
    N        = 8,          # increase (e.g., to 20) for larger networks
    params   = params,
    params_A = params_A,
    params_L = params_L,
    n_iter   = 5000,
    burn_in  = 1000,
    n_alloc  = 50,
    n_cores  = NULL,
    verbose  = TRUE,
    seed     = 42
  )

  # Optional: save output
  # saveRDS(results_s1, file = "results_setting1.rds")
  # write.csv(results_s1$summary, "results_setting1_summary.csv", row.names = FALSE)
}

