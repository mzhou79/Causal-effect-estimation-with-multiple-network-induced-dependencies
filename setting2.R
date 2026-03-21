library(parallel)
library(doParallel)
library(foreach)

# -----------------------------------------------------------------------------
# Group-level network structures (size = 4)
# -----------------------------------------------------------------------------

create_group_network_a_s2 <- function() {
  G <- matrix(0, 4, 4)
  # chain 1-2-3-4 for both interference and dependence
  for (k in 1:3) { G[k, k+1] <- G[k+1, k] <- 1 }
  list(G_interference = G, G_dependence = G, group_size = 4)
}

create_group_network_b_s2 <- function() {
  G <- matrix(0, 4, 4)
  # pairs (1-2), (3-4) for both interference and dependence
  G[1,2] <- G[2,1] <- 1
  G[3,4] <- G[4,3] <- 1
  list(G_interference = G, G_dependence = G, group_size = 4)
}

create_group_network_c_s2 <- function() {
  G_int <- matrix(0, 4, 4)
  G_dep <- matrix(0, 4, 4)
  # interference: pairs only
  G_int[1,2] <- G_int[2,1] <- 1
  G_int[3,4] <- G_int[4,3] <- 1
  # dependence: full chain
  for (k in 1:3) { G_dep[k, k+1] <- G_dep[k+1, k] <- 1 }
  list(G_interference = G_int, G_dependence = G_dep, group_size = 4)
}

# -----------------------------------------------------------------------------
# Weighted energy functions (matching Setting 1)
# Weights are 1/n_neighbors, consistent with Setting 1's energy functions.
# This ensures the true estimands are identical between settings.
# -----------------------------------------------------------------------------

calculate_energy_L_s2 <- function(L_vec, params_L, network) {
  n      <- length(L_vec)
  energy <- sum(L_vec) * params_L$alpha
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (network$G_dependence[i, j] == 1) {
        # weighted: average of the two units' neighbor counts
        n_i   <- sum(network$G_dependence[i, ])
        n_j   <- sum(network$G_dependence[j, ])
        w_ij  <- 0.5 * (1/n_i + 1/n_j)
        energy <- energy + params_L$omega * w_ij * L_vec[i] * L_vec[j]
      }
    }
  }
  energy
}

calculate_energy_A_s2 <- function(A_vec, L_vec, params_A, network) {
  n      <- length(A_vec)
  energy <- 0
  for (i in 1:n) {
    nb_int <- which(network$G_interference[i, ] == 1)
    lp     <- params_A$gamma0 + params_A$gamma1 * L_vec[i]
    if (length(nb_int) > 0)
      lp <- lp + params_A$gamma2 * mean(L_vec[nb_int])   # mean = weighted by 1/n
    energy <- energy + A_vec[i] * lp
  }
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (network$G_dependence[i, j] == 1) {
        n_i  <- sum(network$G_dependence[i, ])
        n_j  <- sum(network$G_dependence[j, ])
        w_ij <- 0.5 * (1/n_i + 1/n_j)
        energy <- energy + params_A$psi * w_ij * A_vec[i] * A_vec[j]
      }
    }
  }
  energy
}

calculate_energy_Y_s2 <- function(Y_vec, A_vec, L_vec, params, network) {
  n      <- length(Y_vec)
  energy <- 0
  for (i in 1:n) {
    nb_int <- which(network$G_interference[i, ] == 1)
    lp     <- params$beta0 + params$beta1 * A_vec[i] + params$beta2 * L_vec[i]
    if (length(nb_int) > 0) {
      lp <- lp + params$beta3 * mean(A_vec[nb_int]) +
                 params$beta4 * mean(L_vec[nb_int])
    }
    energy <- energy + Y_vec[i] * lp
  }
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (network$G_dependence[i, j] == 1) {
        n_i  <- sum(network$G_dependence[i, ])
        n_j  <- sum(network$G_dependence[j, ])
        w_ij <- 0.5 * (1/n_i + 1/n_j)
        energy <- energy + params$theta * w_ij * Y_vec[i] * Y_vec[j]
      }
    }
  }
  energy
}

# -----------------------------------------------------------------------------
# SECTION S2-3: Analytical truth under true model (c)
# Uses the same weighted energy functions and same parameters as Setting 1
# -----------------------------------------------------------------------------

calculate_true_estimands_s2 <- function(params, params_A, params_L) {

  network_true <- create_group_network_c_s2()
  gs           <- 4
  L_cfg        <- generate_binary_configs(gs)
  A_cfg        <- generate_binary_configs(gs)
  Y_cfg        <- generate_binary_configs(gs)

  # f(L)
  fL_unnorm <- apply(L_cfg, 1, function(l)
    exp(calculate_energy_L_s2(l, params_L, network_true)))
  fL <- fL_unnorm / sum(fL_unnorm)

  # pi(a) = sum_l f(a|l) f(l)
  pi_A <- rep(0, nrow(A_cfg))
  for (li in 1:nrow(L_cfg)) {
    fA_unnorm <- apply(A_cfg, 1, function(a)
      exp(calculate_energy_A_s2(a, L_cfg[li,], params_A, network_true)))
    pi_A <- pi_A + fL[li] * fA_unnorm / sum(fA_unnorm)
  }

  # psi_i(a) = sum_l E[Yi | do(A=a), L=l] * f(l)
  calc_psi <- function(a_vec) {
    psi <- rep(0, gs)
    for (li in 1:nrow(L_cfg)) {
      fY_unnorm <- apply(Y_cfg, 1, function(y)
        exp(calculate_energy_Y_s2(y, a_vec, L_cfg[li,], params, network_true)))
      fY    <- fY_unnorm / sum(fY_unnorm)
      EY    <- colSums(Y_cfg * fY)
      psi   <- psi + fL[li] * EY
    }
    psi
  }

  psi_1 <- psi_0 <- rep(0, gs)
  for (i in 1:gs) {
    idx1 <- which(A_cfg[, i] == 1); w1 <- pi_A[idx1] / sum(pi_A[idx1])
    idx0 <- which(A_cfg[, i] == 0); w0 <- pi_A[idx0] / sum(pi_A[idx0])
    tmp1 <- tmp0 <- rep(0, gs)
    for (k in seq_along(idx1)) tmp1 <- tmp1 + w1[k] * calc_psi(A_cfg[idx1[k],])
    for (k in seq_along(idx0)) tmp0 <- tmp0 + w0[k] * calc_psi(A_cfg[idx0[k],])
    psi_1[i] <- tmp1[i]
    psi_0[i] <- tmp0[i]
  }
  psi_00 <- calc_psi(rep(0, gs))

  list(
    DE      = mean(psi_1 - psi_0),
    IE      = mean(psi_0 - psi_00),
    ATE     = mean(psi_1 - psi_00),
    DE_unit = psi_1 - psi_0,
    IE_unit = psi_0 - psi_00,
    psi_1   = psi_1,
    psi_0   = psi_0,
    psi_00  = psi_00
  )
}

# -----------------------------------------------------------------------------
# Gibbs samplers at group level
# -----------------------------------------------------------------------------

# Observational sampler under true graph (c) — used to draw realistic A allocations
gibbs_obs_group_c_s2 <- function(params, params_A, params_L,
                                  n_iter = 4000, burn_in = 1000) {
  net <- create_group_network_c_s2()
  n   <- 4
  L   <- rbinom(n, 1, 0.5); A <- rbinom(n, 1, 0.5); Y <- rbinom(n, 1, 0.5)
  n_keep <- n_iter - burn_in
  Ls <- matrix(0, n_keep, n); As <- matrix(0, n_keep, n); Ys <- matrix(0, n_keep, n)
  ki <- 1
  for (m in 1:n_iter) {
    for (i in sample(1:n)) {
      L1 <- L; L1[i] <- 1; L0 <- L; L0[i] <- 0
      e1 <- calculate_energy_L_s2(L1, params_L, net)
      e0 <- calculate_energy_L_s2(L0, params_L, net)
      me <- max(e0, e1)
      L[i] <- rbinom(1, 1, exp(e1-me) / (exp(e0-me) + exp(e1-me)))

      A1 <- A; A1[i] <- 1; A0 <- A; A0[i] <- 0
      e1 <- calculate_energy_A_s2(A1, L, params_A, net)
      e0 <- calculate_energy_A_s2(A0, L, params_A, net)
      me <- max(e0, e1)
      A[i] <- rbinom(1, 1, exp(e1-me) / (exp(e0-me) + exp(e1-me)))

      Y1 <- Y; Y1[i] <- 1; Y0 <- Y; Y0[i] <- 0
      e1 <- calculate_energy_Y_s2(Y1, A, L, params, net)
      e0 <- calculate_energy_Y_s2(Y0, A, L, params, net)
      me <- max(e0, e1)
      Y[i] <- rbinom(1, 1, exp(e1-me) / (exp(e0-me) + exp(e1-me)))
    }
    if (m > burn_in) { Ls[ki,] <- L; As[ki,] <- A; Ys[ki,] <- Y; ki <- ki+1 }
  }
  list(L_samples = Ls, A_samples = As, Y_samples = Ys)
}

# Interventional sampler under a fitted graph — A is fixed at a_vec
gibbs_intervene_group_s2 <- function(a_vec, params, params_L, network_fit,
                                      n_iter = 4000, burn_in = 1000) {
  n   <- 4
  L   <- rbinom(n, 1, 0.5); Y <- rbinom(n, 1, 0.5)
  n_keep <- n_iter - burn_in
  Ys  <- matrix(0, n_keep, n); ki <- 1
  for (m in 1:n_iter) {
    for (i in sample(1:n)) {
      L1 <- L; L1[i] <- 1; L0 <- L; L0[i] <- 0
      e1 <- calculate_energy_L_s2(L1, params_L, network_fit)
      e0 <- calculate_energy_L_s2(L0, params_L, network_fit)
      me <- max(e0, e1)
      L[i] <- rbinom(1, 1, exp(e1-me) / (exp(e0-me) + exp(e1-me)))

      Y1 <- Y; Y1[i] <- 1; Y0 <- Y; Y0[i] <- 0
      e1 <- calculate_energy_Y_s2(Y1, a_vec, L, params, network_fit)
      e0 <- calculate_energy_Y_s2(Y0, a_vec, L, params, network_fit)
      me <- max(e0, e1)
      Y[i] <- rbinom(1, 1, exp(e1-me) / (exp(e0-me) + exp(e1-me)))
    }
    if (m > burn_in) { Ys[ki,] <- Y; kai <- ki+1 }
  }
  colMeans(Ys)
}

# -----------------------------------------------------------------------------
# Single simulation replicate
# -----------------------------------------------------------------------------

run_single_replicate_s2 <- function(sim_id, network_fit,
                                     params, params_A, params_L,
                                     n_groups, n_alloc,
                                     n_iter_intervene, burn_in_intervene,
                                     seed = NULL) {
  if (!is.null(seed)) set.seed(seed + sim_id)

  # Draw observed A allocations from the TRUE graph (c)
  A_obs <- array(0, dim = c(n_alloc, n_groups, 4))
  for (g in 1:n_groups) {
    obs  <- gibbs_obs_group_c_s2(params, params_A, params_L,
                                  n_iter  = n_alloc * 60,
                                  burn_in = 1000)
    pick <- round(seq(1, nrow(obs$A_samples), length.out = n_alloc))
    A_obs[, g, ] <- obs$A_samples[pick, ]
  }

  # Estimate psi_i(1) and psi_i(0) by averaging over observed allocations,
  # but computing E[Y | do(A = a_obs)] under the FITTED graph
  N       <- n_groups * 4
  psi_1   <- psi_0 <- cnt1 <- cnt0 <- numeric(N)

  for (k in 1:n_alloc) {
    for (g in 1:n_groups) {
      a_vec      <- A_obs[k, g, ]
      EY         <- gibbs_intervene_group_s2(a_vec, params, params_L, network_fit,
                                              n_iter_intervene, burn_in_intervene)
      idx_global <- ((g - 1) * 4 + 1):(g * 4)
      for (i in 1:4) {
        if (a_vec[i] == 1) {
          psi_1[idx_global[i]] <- psi_1[idx_global[i]] + EY[i]
          cnt1[idx_global[i]]  <- cnt1[idx_global[i]]  + 1
        } else {
          psi_0[idx_global[i]] <- psi_0[idx_global[i]] + EY[i]
          cnt0[idx_global[i]]  <- cnt0[idx_global[i]]  + 1
        }
      }
    }
  }

  psi_1 <- psi_1 / pmax(cnt1, 1)
  psi_0 <- psi_0 / pmax(cnt0, 1)

  # Estimate psi_i(0,...,0) under the fitted graph
  psi_00 <- numeric(N)
  for (g in 1:n_groups) {
    EY0        <- gibbs_intervene_group_s2(rep(0, 4), params, params_L, network_fit,
                                            n_iter_intervene, burn_in_intervene)
    idx_global <- ((g - 1) * 4 + 1):(g * 4)
    psi_00[idx_global] <- EY0
  }

  DE  <- mean(psi_1 - psi_0)
  IE  <- mean(psi_0 - psi_00)
  ATE <- mean(psi_1 - psi_00)
  return(list(DE = DE, IE = IE, ATE = ATE))
}

# -----------------------------------------------------------------------------
# Main parallel simulation — Setting 2
# -----------------------------------------------------------------------------

run_simulation_setting2_parallel <- function(
    n_sim    = 500,
    n_groups = 5,
    params   = list(beta0=-1, beta1=0.5, beta2=0.3, beta3=0.4, beta4=0.2, theta=0.3),
    params_A = list(gamma0=0, gamma1=0.3, gamma2=0.2, psi=0.4),
    params_L = list(alpha=-0.5, omega=0.4),
    n_alloc           = 50,
    n_iter_intervene  = 5000,
    burn_in_intervene = 1000,
    n_cores  = NULL,
    verbose  = TRUE,
    seed     = 42) {

  if (is.null(n_cores)) n_cores <- max(1, detectCores() - 1)

  net_a <- create_group_network_a_s2()
  net_b <- create_group_network_b_s2()
  net_c <- create_group_network_c_s2()

  if (verbose) {
    cat("==============================================\n")
    cat("SETTING 2: NETWORK MISSPECIFICATION STUDY\n")
    cat("True DGP: Graph (c) | Fitted: (a), (b), (c)\n")
    cat("n_sim =", n_sim, "| n_cores =", n_cores, "\n")
    cat("==============================================\n\n")
  }

  if (verbose) cat("Calculating analytical truth...\n")
  truth <- calculate_true_estimands_s2(params, params_A, params_L)
  if (verbose) {
    cat("True DE:", round(truth$DE, 4),
        "| True IE:", round(truth$IE, 4),
        "| True ATE:", round(truth$ATE, 4), "\n\n")
  }

  # Functions to export to parallel workers
  fns_to_export <- c(
    "create_group_network_a_s2", "create_group_network_b_s2",
    "create_group_network_c_s2", "generate_binary_configs",
    "calculate_energy_L_s2", "calculate_energy_A_s2", "calculate_energy_Y_s2",
    "gibbs_obs_group_c_s2", "gibbs_intervene_group_s2",
    "run_single_replicate_s2"
  )

  # Run each fitted graph in parallel
  run_parallel_for_graph <- function(network_fit, graph_label) {
    if (verbose) cat("Running graph", graph_label, "in parallel...\n")
    cl <- makeCluster(n_cores); registerDoParallel(cl)
    clusterExport(cl, fns_to_export, envir = environment())
    t0 <- Sys.time()
    res <- foreach(sim = 1:n_sim, .combine = rbind) %dopar% {
      r <- run_single_replicate_s2(
        sim_id            = sim,
        network_fit       = network_fit,
        params            = params,
        params_A          = params_A,
        params_L          = params_L,
        n_groups          = n_groups,
        n_alloc           = n_alloc,
        n_iter_intervene  = n_iter_intervene,
        burn_in_intervene = burn_in_intervene,
        seed              = seed)
      c(DE = r$DE, IE = r$IE, ATE = r$ATE)
    }
    stopCluster(cl)
    elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
    if (verbose) cat("  Graph", graph_label, "done in", elapsed, "minutes\n")
    res
  }

  res_a <- run_parallel_for_graph(net_a, "(a)")
  res_b <- run_parallel_for_graph(net_b, "(b)")
  res_c <- run_parallel_for_graph(net_c, "(c)")

  # Coverage: fraction of estimates within 1.96 SDs of the true value
  coverage <- function(est, tv) mean(abs(est - tv) <= 1.96 * sd(est))

  # Build summary table
  make_rows <- function(res_mat, model_label) {
    do.call(rbind, lapply(c("DE", "IE", "ATE"), function(nm) {
      tv  <- truth[[nm]]
      est <- res_mat[, nm]
      data.frame(
        model    = model_label,
        estimand = nm,
        truth    = round(tv, 4),
        mean     = round(mean(est), 4),
        bias     = round(mean(est) - tv, 4),
        se       = round(sd(est), 4),
        rmse     = round(sqrt(mean((est - tv)^2)), 4),
        coverage = round(coverage(est, tv) * 100, 1),
        stringsAsFactors = FALSE
      )
    }))
  }

  summary_table <- rbind(
    make_rows(res_a, "Graph (a)"),
    make_rows(res_b, "Graph (b)"),
    make_rows(res_c, "Graph (c)")
  )

  if (verbose) {
    cat("\n==============================================\n")
    cat("RESULTS — SETTING 2\n")
    cat("==============================================\n")
    print(summary_table, row.names = FALSE)
  }

  list(
    truth         = truth,
    raw_a         = res_a,
    raw_b         = res_b,
    raw_c         = res_c,
    summary       = summary_table
  )
}

# =============================================================================
# RUN SETTING 2
# =============================================================================

# Same parameters as Setting 1
params   <- list(beta0=-1, beta1=0.5, beta2=0.3, beta3=0.4, beta4=0.2, theta=0.3)
params_A <- list(gamma0=0, gamma1=0.3, gamma2=0.2, psi=0.4)
params_L <- list(alpha=-0.5, omega=0.4)

results_s2 <- run_simulation_setting2_parallel(
  n_sim             = 500,
  n_groups          = 5,
  params            = params,
  params_A          = params_A,
  params_L          = params_L,
  n_alloc           = 50,
  n_iter_intervene  = 5000,
  burn_in_intervene = 1000,
  n_cores           = NULL,
  verbose           = TRUE,
  seed              = 42
)

saveRDS(results_s2, "~/results_s2.rds")
cat("Setting 2 results saved.\n")
