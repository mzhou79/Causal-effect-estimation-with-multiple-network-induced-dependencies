library(parallel)
library(doParallel)
library(foreach)

# =============================================================================
# GRAPH (a)
# =============================================================================

# -----------------------------------------------------------------------------
# Network structure
# -----------------------------------------------------------------------------

create_network_structure_a <- function(N = 20) {
  n_groups   <- N / 4
  group_size <- 4
  G_interference <- matrix(0, N, N)
  G_dependence   <- matrix(0, N, N)
  for (g in 1:n_groups) {
    idx <- ((g - 1) * group_size + 1):(g * group_size)
    for (k in 1:(group_size - 1)) {
      i <- idx[k]; j <- idx[k + 1]
      G_interference[i, j] <- 1; G_interference[j, i] <- 1
      G_dependence[i, j]   <- 1; G_dependence[j, i]   <- 1
    }
  }
  return(list(G_interference = G_interference, G_dependence = G_dependence,
              N = N, n_groups = n_groups, group_size = group_size))
}

# -----------------------------------------------------------------------------
# Helper functions (possible configurations)
# -----------------------------------------------------------------------------

generate_binary_configs <- function(n) {
  configs <- matrix(0, 2^n, n)
  for (j in 1:n)
    configs[, j] <- rep(c(0, 1), each = 2^(n - j), times = 2^(j - 1))
  return(configs)
}

get_chain_neighbors_local <- function(i, group_size = 4) {
  neighbors <- c()
  if (i > 1)          neighbors <- c(neighbors, i - 1)
  if (i < group_size) neighbors <- c(neighbors, i + 1)
  return(neighbors)
}

# -----------------------------------------------------------------------------
# Energy functions
# -----------------------------------------------------------------------------

calculate_group_energy_Y <- function(Y_vec, A_vec, L_vec, params, group_size = 4) {
  energy <- 0
  for (i in 1:group_size) {
    nb      <- get_chain_neighbors_local(i, group_size)
    n_nb    <- length(nb)
    lp      <- params$beta0 + params$beta1 * A_vec[i] + params$beta2 * L_vec[i]
    if (n_nb > 0) {
      w  <- 1 / n_nb
      lp <- lp + params$beta3 * sum(w * A_vec[nb]) + params$beta4 * sum(w * L_vec[nb])
    }
    energy <- energy + Y_vec[i] * lp
  }
  for (i in 1:(group_size - 1)) {
    j     <- i + 1
    w_avg <- 0.5 * (1 / length(get_chain_neighbors_local(i, group_size)) +
                      1 / length(get_chain_neighbors_local(j, group_size)))
    energy <- energy + Y_vec[i] * Y_vec[j] * w_avg * params$theta
  }
  return(energy)
}

calculate_group_energy_A <- function(A_vec, L_vec, params_A, group_size = 4) {
  energy <- 0
  for (i in 1:group_size) {
    nb   <- get_chain_neighbors_local(i, group_size)
    n_nb <- length(nb)
    energy <- energy + A_vec[i] * (params_A$gamma0 + params_A$gamma1 * L_vec[i])
    if (n_nb > 0)
      energy <- energy + A_vec[i] * params_A$gamma2 * sum((1 / n_nb) * L_vec[nb])
  }
  for (i in 1:(group_size - 1)) {
    j     <- i + 1
    w_avg <- 0.5 * (1 / length(get_chain_neighbors_local(i, group_size)) +
                      1 / length(get_chain_neighbors_local(j, group_size)))
    energy <- energy + A_vec[i] * A_vec[j] * w_avg * params_A$psi
  }
  return(energy)
}

calculate_group_energy_L <- function(L_vec, params_L, group_size = 4) {
  energy <- sum(L_vec) * params_L$alpha
  for (i in 1:(group_size - 1)) {
    j     <- i + 1
    w_avg <- 0.5 * (1 / length(get_chain_neighbors_local(i, group_size)) +
                      1 / length(get_chain_neighbors_local(j, group_size)))
    energy <- energy + L_vec[i] * L_vec[j] * w_avg * params_L$omega
  }
  return(energy)
}

# -----------------------------------------------------------------------------
# Conditional probability functions
# -----------------------------------------------------------------------------

prob_Y_given_rest <- function(i, Y, A, L, params, network) {
  gs      <- network$group_size
  g       <- ceiling(i / gs)
  i_loc   <- i - (g - 1) * gs
  idx     <- ((g - 1) * gs + 1):(g * gs)
  Yg <- Y[idx]; Ag <- A[idx]; Lg <- L[idx]
  Y1 <- Yg; Y1[i_loc] <- 1
  Y0 <- Yg; Y0[i_loc] <- 0
  e1 <- calculate_group_energy_Y(Y1, Ag, Lg, params, gs)
  e0 <- calculate_group_energy_Y(Y0, Ag, Lg, params, gs)
  me <- max(e0, e1)
  return(exp(e1 - me) / (exp(e0 - me) + exp(e1 - me)))
}

prob_A_given_rest <- function(i, A, L, params_A, network) {
  gs    <- network$group_size
  g     <- ceiling(i / gs)
  i_loc <- i - (g - 1) * gs
  idx   <- ((g - 1) * gs + 1):(g * gs)
  Ag <- A[idx]; Lg <- L[idx]
  A1 <- Ag; A1[i_loc] <- 1
  A0 <- Ag; A0[i_loc] <- 0
  e1 <- calculate_group_energy_A(A1, Lg, params_A, gs)
  e0 <- calculate_group_energy_A(A0, Lg, params_A, gs)
  me <- max(e0, e1)
  return(exp(e1 - me) / (exp(e0 - me) + exp(e1 - me)))
}

prob_L_given_rest <- function(i, L, params_L, network) {
  gs    <- network$group_size
  g     <- ceiling(i / gs)
  i_loc <- i - (g - 1) * gs
  idx   <- ((g - 1) * gs + 1):(g * gs)
  Lg <- L[idx]
  L1 <- Lg; L1[i_loc] <- 1
  L0 <- Lg; L0[i_loc] <- 0
  e1 <- calculate_group_energy_L(L1, params_L, gs)
  e0 <- calculate_group_energy_L(L0, params_L, gs)
  me <- max(e0, e1)
  return(exp(e1 - me) / (exp(e0 - me) + exp(e1 - me)))
}

# -----------------------------------------------------------------------------
# Analytical estimands
# ATE(alpha) = DE(alpha) + IE(alpha) = psi_1 - psi_0_0
# -----------------------------------------------------------------------------

calculate_psi_analytical_group <- function(i, a_vec, params, params_L, group_size = 4) {
  L_cfg <- generate_binary_configs(group_size)
  Y_cfg <- generate_binary_configs(group_size)
  fL    <- exp(apply(L_cfg, 1, function(l) calculate_group_energy_L(l, params_L, group_size)))
  fL    <- fL / sum(fL)
  EYi   <- sapply(1:nrow(L_cfg), function(li) {
    fY <- exp(apply(Y_cfg, 1, function(y) calculate_group_energy_Y(y, a_vec, L_cfg[li, ], params, group_size)))
    sum(Y_cfg[, i] * fY / sum(fY))
  })
  return(sum(EYi * fL))
}

calculate_prob_A_given_L <- function(A_vec, L_vec, params_A, group_size = 4) {
  A_cfg  <- generate_binary_configs(group_size)
  fA     <- exp(apply(A_cfg, 1, function(a) calculate_group_energy_A(a, L_vec, params_A, group_size)))
  a_idx  <- which(apply(A_cfg, 1, function(x) all(x == A_vec)))
  return(fA[a_idx] / sum(fA))
}

calculate_marginal_pi_A <- function(A_vec, params_A, params_L, group_size = 4) {
  L_cfg <- generate_binary_configs(group_size)
  fL    <- exp(apply(L_cfg, 1, function(l) calculate_group_energy_L(l, params_L, group_size)))
  fL    <- fL / sum(fL)
  sum(sapply(1:nrow(L_cfg), function(li)
    calculate_prob_A_given_L(A_vec, L_cfg[li, ], params_A, group_size) * fL[li]))
}

calculate_estimands_analytical <- function(network, params, params_A, params_L) {
  gs       <- network$group_size
  n_groups <- network$n_groups
  A_cfg    <- generate_binary_configs(gs)
  cat("  Computing allocation probabilities pi(a)...\n")
  pi_A <- sapply(1:nrow(A_cfg), function(k)
    calculate_marginal_pi_A(A_cfg[k, ], params_A, params_L, gs))
  cat("  Allocation probabilities sum to:", round(sum(pi_A), 6), "\n")
  cat("  Computing analytical psi for each position...\n")
  psi_1 <- psi_0 <- psi_0_0 <- numeric(gs)
  for (i in 1:gs) {
    cat("    Position", i, "\n")
    cfg1    <- which(A_cfg[, i] == 1); pi1 <- pi_A[cfg1] / sum(pi_A[cfg1])
    psi_1[i] <- sum(pi1 * sapply(cfg1, function(k)
      calculate_psi_analytical_group(i, A_cfg[k, ], params, params_L, gs)))
    cfg0    <- which(A_cfg[, i] == 0); pi0 <- pi_A[cfg0] / sum(pi_A[cfg0])
    psi_0[i] <- sum(pi0 * sapply(cfg0, function(k)
      calculate_psi_analytical_group(i, A_cfg[k, ], params, params_L, gs)))
    psi_0_0[i] <- calculate_psi_analytical_group(i, rep(0, gs), params, params_L, gs)
  }
  DE  <- mean(psi_1 - psi_0)
  IE  <- mean(psi_0 - psi_0_0)
  ATE <- mean(psi_1 - psi_0_0)   # ATE = DE + IE
  return(list(DE = DE, IE = IE, ATE = ATE,
              psi_1 = rep(psi_1, n_groups), psi_0 = rep(psi_0, n_groups),
              psi_0_0 = rep(psi_0_0, n_groups),
              psi_1_group = psi_1, psi_0_group = psi_0, psi_0_0_group = psi_0_0,
              pi_A = pi_A))
}

# -----------------------------------------------------------------------------
# Gibbs samplers
# -----------------------------------------------------------------------------

gibbs_sampler_observational <- function(network, params, params_A, params_L,
                                        n_iter = 5000, burn_in = 1000) {
  N <- network$N
  L <- rbinom(N, 1, 0.5); A <- rbinom(N, 1, 0.5); Y <- rbinom(N, 1, 0.5)
  n_s <- n_iter - burn_in
  Ls <- matrix(0, n_s, N); As <- matrix(0, n_s, N); Ys <- matrix(0, n_s, N)
  si <- 1
  for (m in 1:n_iter) {
    for (i in sample(1:N)) {
      L[i] <- rbinom(1, 1, prob_L_given_rest(i, L, params_L, network))
      A[i] <- rbinom(1, 1, prob_A_given_rest(i, A, L, params_A, network))
      Y[i] <- rbinom(1, 1, prob_Y_given_rest(i, Y, A, L, params, network))
    }
    if (m > burn_in) { Ls[si, ] <- L; As[si, ] <- A; Ys[si, ] <- Y; si <- si + 1 }
  }
  return(list(L_samples = Ls, A_samples = As, Y_samples = Ys))
}

gibbs_sampler_interventional <- function(a_vec, network, params, params_L,
                                         n_iter = 5000, burn_in = 1000) {
  N <- network$N
  L <- rbinom(N, 1, 0.5); Y <- rbinom(N, 1, 0.5)
  n_s <- n_iter - burn_in
  Ls <- matrix(0, n_s, N); Ys <- matrix(0, n_s, N)
  si <- 1
  for (m in 1:n_iter) {
    for (i in sample(1:N)) {
      L[i] <- rbinom(1, 1, prob_L_given_rest(i, L, params_L, network))
      Y[i] <- rbinom(1, 1, prob_Y_given_rest(i, Y, a_vec, L, params, network))
    }
    if (m > burn_in) { Ls[si, ] <- L; Ys[si, ] <- Y; si <- si + 1 }
  }
  return(list(L_samples = Ls, Y_samples = Ys))
}

# -----------------------------------------------------------------------------
# Single simulation run — returns DE, IE, ATE
# -----------------------------------------------------------------------------

run_single_simulation_a <- function(sim_id, network, params, params_A, params_L,
                                    n_iter, burn_in, n_alloc, seed = NULL) {
  if (!is.null(seed)) set.seed(seed + sim_id)
  N <- network$N
  obs <- gibbs_sampler_observational(network, params, params_A, params_L,
                                     n_iter = n_alloc * 50, burn_in = 500)
  aidx <- seq(1, nrow(obs$A_samples), length.out = n_alloc)
  psi_1 <- psi_0 <- numeric(N)
  cnt1  <- cnt0  <- numeric(N)
  for (k in 1:n_alloc) {
    a_vec   <- obs$A_samples[round(aidx[k]), ]
    Ymeans  <- colMeans(gibbs_sampler_interventional(a_vec, network, params, params_L,
                                                     n_iter, burn_in)$Y_samples)
    for (i in 1:N) {
      if (a_vec[i] == 1) { psi_1[i] <- psi_1[i] + Ymeans[i]; cnt1[i] <- cnt1[i] + 1 }
      else               { psi_0[i] <- psi_0[i] + Ymeans[i]; cnt0[i] <- cnt0[i] + 1 }
    }
  }
  psi_1   <- psi_1 / pmax(cnt1, 1)
  psi_0   <- psi_0 / pmax(cnt0, 1)
  psi_0_0 <- colMeans(gibbs_sampler_interventional(rep(0, N), network, params, params_L,
                                                   n_iter, burn_in)$Y_samples)
  DE  <- mean(psi_1 - psi_0)
  IE  <- mean(psi_0 - psi_0_0)
  ATE <- mean(psi_1 - psi_0_0)   # = DE + IE
  return(list(DE = DE, IE = IE, ATE = ATE))
}

# -----------------------------------------------------------------------------
# Main parallel simulation — graph (a)
# -----------------------------------------------------------------------------

run_simulation_a_parallel <- function(
    n_sim = 500, N = 20,
    params   = list(beta0=-1, beta1=0.5, beta2=0.3, beta3=0.4, beta4=0.2, theta=0.3),
    params_A = list(gamma0=0, gamma1=0.3, gamma2=0.2, psi=0.4),
    params_L = list(alpha=-0.5, omega=0.4),
    n_iter = 5000, burn_in = 1000, n_alloc = 50,
    n_cores = NULL, verbose = TRUE, seed = 42) {
  
  if (is.null(n_cores)) n_cores <- max(1, detectCores() - 1)
  network <- create_network_structure_a(N)
  
  if (verbose) {
    cat("===========================================\n")
    cat("PARALLEL Simulation for Chain Graph (a)\n")
    cat("Structure: Chain 1-2-3-4 for BOTH\n")
    cat("           interference and dependence\n")
    cat("N =", N, "| n_sim =", n_sim, "| n_cores =", n_cores, "\n")
    cat("===========================================\n\n")
  }
  
  if (verbose) cat("Calculating analytical truth...\n")
  truth <- calculate_estimands_analytical(network, params, params_A, params_L)
  if (verbose) {
    cat("\nTrue DE:", round(truth$DE, 4),
        "| True IE:", round(truth$IE, 4),
        "| True ATE:", round(truth$ATE, 4), "\n\n")
  }
  
  cl <- makeCluster(n_cores); registerDoParallel(cl)
  clusterExport(cl, c(
    "create_network_structure_a", "generate_binary_configs",
    "get_chain_neighbors_local",
    "calculate_group_energy_Y", "calculate_group_energy_A", "calculate_group_energy_L",
    "prob_Y_given_rest", "prob_A_given_rest", "prob_L_given_rest",
    "gibbs_sampler_observational", "gibbs_sampler_interventional",
    "run_single_simulation_a"), envir = environment())
  
  if (verbose) cat("Running", n_sim, "simulations in parallel...\n")
  t0 <- Sys.time()
  res_list <- foreach(sim = 1:n_sim, .combine = rbind) %dopar% {
    r <- run_single_simulation_a(sim, network, params, params_A, params_L,
                                 n_iter, burn_in, n_alloc, seed)
    c(DE = r$DE, IE = r$IE, ATE = r$ATE)
  }
  stopCluster(cl)
  elapsed <- difftime(Sys.time(), t0, units = "mins")
  if (verbose) cat("Completed in", round(as.numeric(elapsed), 2), "minutes\n\n")
  
  DE_est  <- res_list[, "DE"]
  IE_est  <- res_list[, "IE"]
  ATE_est <- res_list[, "ATE"]
  
  # Coverage: fraction of estimates within 1.96 SDs of the true value
  coverage <- function(est, truth_val)
    mean(abs(est - truth_val) <= 1.96 * sd(est))
  
  results <- list(
    true_DE = truth$DE, true_IE = truth$IE, true_ATE = truth$ATE,
    DE_estimates = DE_est, IE_estimates = IE_est, ATE_estimates = ATE_est,
    DE_mean  = mean(DE_est),  DE_bias  = mean(DE_est)  - truth$DE,
    DE_se    = sd(DE_est),    DE_rmse  = sqrt(mean((DE_est  - truth$DE)^2)),
    IE_mean  = mean(IE_est),  IE_bias  = mean(IE_est)  - truth$IE,
    IE_se    = sd(IE_est),    IE_rmse  = sqrt(mean((IE_est  - truth$IE)^2)),
    ATE_mean = mean(ATE_est), ATE_bias = mean(ATE_est) - truth$ATE,
    ATE_se   = sd(ATE_est),   ATE_rmse = sqrt(mean((ATE_est - truth$ATE)^2)),
    DE_coverage  = coverage(DE_est,  truth$DE),
    IE_coverage  = coverage(IE_est,  truth$IE),
    ATE_coverage = coverage(ATE_est, truth$ATE),
    elapsed_time_mins = as.numeric(elapsed), n_cores_used = n_cores)
  
  if (verbose) {
    cat("===========================================\n")
    cat("RESULTS FOR GRAPH (a)\n")
    cat("===========================================\n")
    for (nm in c("DE", "IE", "ATE")) {
      est  <- results[[paste0(nm, "_estimates")]]
      cat(sprintf("\n%s:\n  True: %.4f | Mean: %.4f | Bias: %.4f | SE: %.4f | RMSE: %.4f | 95%% Cov: %.1f%%\n",
                  nm, results[[paste0("true_", nm)]], results[[paste0(nm, "_mean")]],
                  results[[paste0(nm, "_bias")]], results[[paste0(nm, "_se")]],
                  results[[paste0(nm, "_rmse")]], results[[paste0(nm, "_coverage")]] * 100))
    }
    cat("\nComputation time:", round(results$elapsed_time_mins, 2), "minutes\n")
  }
  return(results)
}


# =============================================================================
# GRAPH (b)
# =============================================================================

# -----------------------------------------------------------------------------
# Network structure
# -----------------------------------------------------------------------------

create_network_structure_b <- function(N = 20) {
  n_groups <- N / 4; group_size <- 4
  G_interference <- matrix(0, N, N); G_dependence <- matrix(0, N, N)
  for (g in 1:n_groups) {
    idx <- ((g - 1) * group_size + 1):(g * group_size)
    for (pair in list(c(1,2), c(3,4))) {
      i <- idx[pair[1]]; j <- idx[pair[2]]
      G_interference[i,j] <- G_interference[j,i] <- 1
      G_dependence[i,j]   <- G_dependence[j,i]   <- 1
    }
  }
  return(list(G_interference = G_interference, G_dependence = G_dependence,
              N = N, n_groups = n_groups, group_size = group_size))
}

get_pair_neighbors_local <- function(i, group_size = 4) {
  if (i == 1) return(2); if (i == 2) return(1)
  if (i == 3) return(4); if (i == 4) return(3)
  return(c())
}

get_interference_neighbors <- function(i, G) which(G[i, ] == 1)
get_dependence_neighbors   <- function(i, G) which(G[i, ] == 1)

# -----------------------------------------------------------------------------
# Energy functions
# -----------------------------------------------------------------------------

calculate_group_energy_Y_b <- function(Y_vec, A_vec, L_vec, params, group_size = 4) {
  energy <- 0
  for (i in 1:group_size) {
    nb <- get_pair_neighbors_local(i, group_size); n_nb <- length(nb)
    lp <- params$beta0 + params$beta1 * A_vec[i] + params$beta2 * L_vec[i]
    if (n_nb > 0) lp <- lp + params$beta3 * A_vec[nb] / n_nb + params$beta4 * L_vec[nb] / n_nb
    energy <- energy + Y_vec[i] * lp
  }
  for (pair in list(c(1,2), c(3,4))) {
    i <- pair[1]; j <- pair[2]
    w <- 0.5 * (1/length(get_pair_neighbors_local(i,group_size)) +
                  1/length(get_pair_neighbors_local(j,group_size)))
    energy <- energy + Y_vec[i] * Y_vec[j] * w * params$theta
  }
  return(energy)
}

calculate_group_energy_A_b <- function(A_vec, L_vec, params_A, group_size = 4) {
  energy <- 0
  for (i in 1:group_size) {
    nb <- get_pair_neighbors_local(i, group_size); n_nb <- length(nb)
    energy <- energy + A_vec[i] * (params_A$gamma0 + params_A$gamma1 * L_vec[i])
    if (n_nb > 0) energy <- energy + A_vec[i] * params_A$gamma2 * L_vec[nb] / n_nb
  }
  for (pair in list(c(1,2), c(3,4))) {
    i <- pair[1]; j <- pair[2]
    w <- 0.5 * (1/length(get_pair_neighbors_local(i,group_size)) +
                  1/length(get_pair_neighbors_local(j,group_size)))
    energy <- energy + A_vec[i] * A_vec[j] * w * params_A$psi
  }
  return(energy)
}

calculate_group_energy_L_b <- function(L_vec, params_L, group_size = 4) {
  energy <- sum(L_vec) * params_L$alpha
  for (pair in list(c(1,2), c(3,4))) {
    i <- pair[1]; j <- pair[2]
    w <- 0.5 * (1/length(get_pair_neighbors_local(i,group_size)) +
                  1/length(get_pair_neighbors_local(j,group_size)))
    energy <- energy + L_vec[i] * L_vec[j] * w * params_L$omega
  }
  return(energy)
}

# -----------------------------------------------------------------------------
# Conditional probability functions
# -----------------------------------------------------------------------------

local_energy_Y_b <- function(i, y_i, Y, A, L, params, network) {
  G_int <- network$G_interference; G_dep <- network$G_dependence
  int_nb <- get_interference_neighbors(i, G_int)
  dep_nb <- get_dependence_neighbors(i, G_dep)
  lp <- params$beta0 + params$beta1 * A[i] + params$beta2 * L[i]
  if (length(int_nb) > 0) lp <- lp +
    params$beta3 * mean(A[int_nb]) + params$beta4 * mean(L[int_nb])
  energy <- y_i * lp
  if (length(dep_nb) > 0)
    energy <- energy + sum(y_i * Y[dep_nb] * (1/length(dep_nb)) * params$theta)
  return(energy)
}

local_energy_A_b <- function(i, a_i, A, L, params_A, network) {
  G_int <- network$G_interference; G_dep <- network$G_dependence
  int_nb <- get_interference_neighbors(i, G_int)
  dep_nb <- get_dependence_neighbors(i, G_dep)
  energy <- a_i * (params_A$gamma0 + params_A$gamma1 * L[i])
  if (length(int_nb) > 0) energy <- energy + a_i * params_A$gamma2 * mean(L[int_nb])
  if (length(dep_nb) > 0) energy <- energy + sum(a_i * A[dep_nb] * (1/length(dep_nb)) * params_A$psi)
  return(energy)
}

local_energy_L_b <- function(i, l_i, L, params_L, network) {
  dep_nb <- get_dependence_neighbors(i, network$G_dependence)
  energy <- l_i * params_L$alpha
  if (length(dep_nb) > 0) energy <- energy + sum(l_i * L[dep_nb] * (1/length(dep_nb)) * params_L$omega)
  return(energy)
}

prob_Y_given_rest_b <- function(i, Y, A, L, params, network) {
  e1 <- local_energy_Y_b(i, 1, replace(Y, i, 1), A, L, params, network)
  e0 <- local_energy_Y_b(i, 0, replace(Y, i, 0), A, L, params, network)
  me <- max(e0, e1); exp(e1-me) / (exp(e0-me) + exp(e1-me))
}

prob_A_given_rest_b <- function(i, A, L, params_A, network) {
  e1 <- local_energy_A_b(i, 1, replace(A, i, 1), L, params_A, network)
  e0 <- local_energy_A_b(i, 0, replace(A, i, 0), L, params_A, network)
  me <- max(e0, e1); exp(e1-me) / (exp(e0-me) + exp(e1-me))
}

prob_L_given_rest_b <- function(i, L, params_L, network) {
  e1 <- local_energy_L_b(i, 1, replace(L, i, 1), params_L, network)
  e0 <- local_energy_L_b(i, 0, replace(L, i, 0), params_L, network)
  me <- max(e0, e1); exp(e1-me) / (exp(e0-me) + exp(e1-me))
}

# -----------------------------------------------------------------------------
# Analytical estimands — graph (b)
# -----------------------------------------------------------------------------

calculate_psi_analytical_group_b <- function(i, a_vec, params, params_L, group_size = 4) {
  L_cfg <- generate_binary_configs(group_size)
  Y_cfg <- generate_binary_configs(group_size)
  fL    <- exp(apply(L_cfg, 1, function(l) calculate_group_energy_L_b(l, params_L, group_size)))
  fL    <- fL / sum(fL)
  EYi   <- sapply(1:nrow(L_cfg), function(li) {
    fY <- exp(apply(Y_cfg, 1, function(y) calculate_group_energy_Y_b(y, a_vec, L_cfg[li,], params, group_size)))
    sum(Y_cfg[, i] * fY / sum(fY))
  })
  return(sum(EYi * fL))
}

calculate_prob_A_given_L_b <- function(A_vec, L_vec, params_A, group_size = 4) {
  A_cfg <- generate_binary_configs(group_size)
  fA    <- exp(apply(A_cfg, 1, function(a) calculate_group_energy_A_b(a, L_vec, params_A, group_size)))
  a_idx <- which(apply(A_cfg, 1, function(x) all(x == A_vec)))
  return(fA[a_idx] / sum(fA))
}

calculate_marginal_pi_A_b <- function(A_vec, params_A, params_L, group_size = 4) {
  L_cfg <- generate_binary_configs(group_size)
  fL    <- exp(apply(L_cfg, 1, function(l) calculate_group_energy_L_b(l, params_L, group_size)))
  fL    <- fL / sum(fL)
  sum(sapply(1:nrow(L_cfg), function(li)
    calculate_prob_A_given_L_b(A_vec, L_cfg[li,], params_A, group_size) * fL[li]))
}

calculate_estimands_analytical_b <- function(network, params, params_A, params_L) {
  gs    <- network$group_size; n_groups <- network$n_groups
  A_cfg <- generate_binary_configs(gs)
  cat("  Computing allocation probabilities pi(a)...\n")
  pi_A  <- sapply(1:nrow(A_cfg), function(k)
    calculate_marginal_pi_A_b(A_cfg[k,], params_A, params_L, gs))
  cat("  Allocation probabilities sum to:", round(sum(pi_A), 6), "\n")
  psi_1 <- psi_0 <- psi_0_0 <- numeric(gs)
  cat("  Computing analytical psi for each position...\n")
  for (i in 1:gs) {
    cat("    Position", i, "\n")
    cfg1 <- which(A_cfg[,i]==1); pi1 <- pi_A[cfg1]/sum(pi_A[cfg1])
    psi_1[i] <- sum(pi1 * sapply(cfg1, function(k)
      calculate_psi_analytical_group_b(i, A_cfg[k,], params, params_L, gs)))
    cfg0 <- which(A_cfg[,i]==0); pi0 <- pi_A[cfg0]/sum(pi_A[cfg0])
    psi_0[i] <- sum(pi0 * sapply(cfg0, function(k)
      calculate_psi_analytical_group_b(i, A_cfg[k,], params, params_L, gs)))
    psi_0_0[i] <- calculate_psi_analytical_group_b(i, rep(0,gs), params, params_L, gs)
  }
  DE  <- mean(psi_1 - psi_0)
  IE  <- mean(psi_0 - psi_0_0)
  ATE <- mean(psi_1 - psi_0_0)   # ATE = DE + IE
  return(list(DE=DE, IE=IE, ATE=ATE,
              psi_1=rep(psi_1,n_groups), psi_0=rep(psi_0,n_groups),
              psi_0_0=rep(psi_0_0,n_groups),
              psi_1_group=psi_1, psi_0_group=psi_0, psi_0_0_group=psi_0_0, pi_A=pi_A))
}

# -----------------------------------------------------------------------------
# Gibbs samplers — graph (b)
# -----------------------------------------------------------------------------

gibbs_sampler_observational_b <- function(network, params, params_A, params_L,
                                          n_iter=5000, burn_in=1000) {
  N <- network$N
  L <- rbinom(N,1,0.5); A <- rbinom(N,1,0.5); Y <- rbinom(N,1,0.5)
  n_s <- n_iter - burn_in
  Ls <- matrix(0,n_s,N); As <- matrix(0,n_s,N); Ys <- matrix(0,n_s,N); si <- 1
  for (m in 1:n_iter) {
    for (i in sample(1:N)) {
      L[i] <- rbinom(1,1,prob_L_given_rest_b(i,L,params_L,network))
      A[i] <- rbinom(1,1,prob_A_given_rest_b(i,A,L,params_A,network))
      Y[i] <- rbinom(1,1,prob_Y_given_rest_b(i,Y,A,L,params,network))
    }
    if (m > burn_in) { Ls[si,]<-L; As[si,]<-A; Ys[si,]<-Y; si<-si+1 }
  }
  return(list(L_samples=Ls, A_samples=As, Y_samples=Ys))
}

gibbs_sampler_interventional_b <- function(a_vec, network, params, params_L,
                                           n_iter=5000, burn_in=1000) {
  N <- network$N
  L <- rbinom(N,1,0.5); Y <- rbinom(N,1,0.5)
  n_s <- n_iter - burn_in
  Ls <- matrix(0,n_s,N); Ys <- matrix(0,n_s,N); si <- 1
  for (m in 1:n_iter) {
    for (i in sample(1:N)) {
      L[i] <- rbinom(1,1,prob_L_given_rest_b(i,L,params_L,network))
      Y[i] <- rbinom(1,1,prob_Y_given_rest_b(i,Y,a_vec,L,params,network))
    }
    if (m > burn_in) { Ls[si,]<-L; Ys[si,]<-Y; si<-si+1 }
  }
  return(list(L_samples=Ls, Y_samples=Ys))
}

# -----------------------------------------------------------------------------
# Single simulation run — graph (b), returns DE, IE, ATE
# -----------------------------------------------------------------------------

run_single_simulation_b <- function(sim_id, network, params, params_A, params_L,
                                    n_iter, burn_in, n_alloc, seed=NULL) {
  if (!is.null(seed)) set.seed(seed + sim_id)
  N   <- network$N
  obs <- gibbs_sampler_observational_b(network, params, params_A, params_L,
                                       n_iter = n_alloc * 50, burn_in = 500)
  aidx  <- seq(1, nrow(obs$A_samples), length.out = n_alloc)
  psi_1 <- psi_0 <- cnt1 <- cnt0 <- numeric(N)
  for (k in 1:n_alloc) {
    a_vec  <- obs$A_samples[round(aidx[k]), ]
    Ymeans <- colMeans(gibbs_sampler_interventional_b(a_vec, network, params, params_L,
                                                      n_iter, burn_in)$Y_samples)
    for (i in 1:N) {
      if (a_vec[i]==1) { psi_1[i]<-psi_1[i]+Ymeans[i]; cnt1[i]<-cnt1[i]+1 }
      else             { psi_0[i]<-psi_0[i]+Ymeans[i]; cnt0[i]<-cnt0[i]+1 }
    }
  }
  psi_1   <- psi_1 / pmax(cnt1, 1)
  psi_0   <- psi_0 / pmax(cnt0, 1)
  psi_0_0 <- colMeans(gibbs_sampler_interventional_b(rep(0,N), network, params, params_L,
                                                     n_iter, burn_in)$Y_samples)
  DE  <- mean(psi_1 - psi_0)
  IE  <- mean(psi_0 - psi_0_0)
  ATE <- mean(psi_1 - psi_0_0)   # ATE = DE + IE
  return(list(DE=DE, IE=IE, ATE=ATE))
}

# -----------------------------------------------------------------------------
# Main parallel simulation — graph (b)
# -----------------------------------------------------------------------------

run_simulation_b_parallel <- function(
    n_sim=500, N=20,
    params   = list(beta0=-1, beta1=0.5, beta2=0.3, beta3=0.4, beta4=0.2, theta=0.3),
    params_A = list(gamma0=0, gamma1=0.3, gamma2=0.2, psi=0.4),
    params_L = list(alpha=-0.5, omega=0.4),
    n_iter=5000, burn_in=1000, n_alloc=50,
    n_cores=NULL, verbose=TRUE, seed=42) {
  
  if (is.null(n_cores)) n_cores <- max(1, detectCores() - 1)
  network <- create_network_structure_b(N)
  
  if (verbose) {
    cat("===========================================\n")
    cat("PARALLEL Simulation for Chain Graph (b)\n")
    cat("Structure: Disconnected pairs (1-2),(3-4)\n")
    cat("N =", N, "| n_sim =", n_sim, "| n_cores =", n_cores, "\n")
    cat("===========================================\n\n")
  }
  
  if (verbose) cat("Calculating analytical truth...\n")
  truth <- calculate_estimands_analytical_b(network, params, params_A, params_L)
  if (verbose) cat("\nTrue DE:", round(truth$DE,4),
                   "| True IE:", round(truth$IE,4),
                   "| True ATE:", round(truth$ATE,4), "\n\n")
  
  cl <- makeCluster(n_cores); registerDoParallel(cl)
  clusterExport(cl, c(
    "create_network_structure_b", "generate_binary_configs",
    "get_pair_neighbors_local", "get_interference_neighbors", "get_dependence_neighbors",
    "calculate_group_energy_Y_b", "calculate_group_energy_A_b", "calculate_group_energy_L_b",
    "local_energy_Y_b", "local_energy_A_b", "local_energy_L_b",
    "prob_Y_given_rest_b", "prob_A_given_rest_b", "prob_L_given_rest_b",
    "gibbs_sampler_observational_b", "gibbs_sampler_interventional_b",
    "run_single_simulation_b"), envir = environment())
  
  if (verbose) cat("Running", n_sim, "simulations in parallel...\n")
  t0 <- Sys.time()
  res_list <- foreach(sim = 1:n_sim, .combine = rbind) %dopar% {
    r <- run_single_simulation_b(sim, network, params, params_A, params_L,
                                 n_iter, burn_in, n_alloc, seed)
    c(DE = r$DE, IE = r$IE, ATE = r$ATE)
  }
  stopCluster(cl)
  elapsed <- difftime(Sys.time(), t0, units = "mins")
  if (verbose) cat("Completed in", round(as.numeric(elapsed), 2), "minutes\n\n")
  
  DE_est  <- res_list[, "DE"]
  IE_est  <- res_list[, "IE"]
  ATE_est <- res_list[, "ATE"]
  
  coverage <- function(est, tv) mean(abs(est - tv) <= 1.96 * sd(est))
  
  results <- list(
    true_DE=truth$DE, true_IE=truth$IE, true_ATE=truth$ATE,
    DE_estimates=DE_est, IE_estimates=IE_est, ATE_estimates=ATE_est,
    DE_mean=mean(DE_est),   DE_bias=mean(DE_est)-truth$DE,
    DE_se=sd(DE_est),       DE_rmse=sqrt(mean((DE_est-truth$DE)^2)),
    IE_mean=mean(IE_est),   IE_bias=mean(IE_est)-truth$IE,
    IE_se=sd(IE_est),       IE_rmse=sqrt(mean((IE_est-truth$IE)^2)),
    ATE_mean=mean(ATE_est), ATE_bias=mean(ATE_est)-truth$ATE,
    ATE_se=sd(ATE_est),     ATE_rmse=sqrt(mean((ATE_est-truth$ATE)^2)),
    DE_coverage=coverage(DE_est,truth$DE),
    IE_coverage=coverage(IE_est,truth$IE),
    ATE_coverage=coverage(ATE_est,truth$ATE),
    elapsed_time_mins=as.numeric(elapsed), n_cores_used=n_cores)
  
  if (verbose) {
    cat("===========================================\n")
    cat("RESULTS FOR GRAPH (b)\n")
    cat("===========================================\n")
    for (nm in c("DE","IE","ATE"))
      cat(sprintf("\n%s:\n  True: %.4f | Mean: %.4f | Bias: %.4f | SE: %.4f | RMSE: %.4f | 95%% Cov: %.1f%%\n",
                  nm, results[[paste0("true_",nm)]], results[[paste0(nm,"_mean")]],
                  results[[paste0(nm,"_bias")]], results[[paste0(nm,"_se")]],
                  results[[paste0(nm,"_rmse")]], results[[paste0(nm,"_coverage")]]*100))
    cat("\nComputation time:", round(results$elapsed_time_mins, 2), "minutes\n")
    cat("Cores used:", results$n_cores_used, "\n")
  }
  return(results)
}


# =============================================================================
# GRAPH (c)
# =============================================================================

# -----------------------------------------------------------------------------
# Network structure
# -----------------------------------------------------------------------------

create_network_structure_c <- function(N = 20) {
  n_groups <- N / 4; group_size <- 4
  G_interference <- matrix(0, N, N); G_dependence <- matrix(0, N, N)
  for (g in 1:n_groups) {
    idx <- ((g - 1) * group_size + 1):(g * group_size)
    for (pair in list(c(1,2), c(3,4))) {
      i <- idx[pair[1]]; j <- idx[pair[2]]
      G_interference[i,j] <- G_interference[j,i] <- 1
    }
    for (k in 1:(group_size-1)) {
      i <- idx[k]; j <- idx[k+1]
      G_dependence[i,j] <- G_dependence[j,i] <- 1
    }
  }
  return(list(G_interference=G_interference, G_dependence=G_dependence,
              N=N, n_groups=n_groups, group_size=group_size))
}

get_interference_neighbors_c <- function(i_local, group_size=4) {
  if (i_local==1) return(2); if (i_local==2) return(1)
  if (i_local==3) return(4); if (i_local==4) return(3)
  return(c())
}

get_dependence_neighbors_c <- function(i_local, group_size=4) {
  nb <- c()
  if (i_local > 1)          nb <- c(nb, i_local-1)
  if (i_local < group_size) nb <- c(nb, i_local+1)
  return(nb)
}

# -----------------------------------------------------------------------------
# Energy functions — graph (c)
# Ni ∩ g(i) = interference neighbors  →  beta3/beta4 and A-A psi terms
# Ni         = dependence neighbors   →  theta Y-Y and L-L omega terms
# -----------------------------------------------------------------------------

calculate_group_energy_Y_c <- function(Y_vec, A_vec, L_vec, params, group_size=4) {
  energy <- 0
  for (i in 1:group_size) {
    int_nb <- get_interference_neighbors_c(i, group_size)
    lp <- params$beta0 + params$beta1*A_vec[i] + params$beta2*L_vec[i]
    if (length(int_nb)>0) lp <- lp +
      params$beta3*mean(A_vec[int_nb]) + params$beta4*mean(L_vec[int_nb])
    energy <- energy + Y_vec[i] * lp
  }
  for (i in 1:(group_size-1)) {
    j <- i+1
    dep_i <- get_dependence_neighbors_c(i, group_size)
    dep_j <- get_dependence_neighbors_c(j, group_size)
    w <- 0.5*(1/length(dep_i) + 1/length(dep_j))
    energy <- energy + Y_vec[i]*Y_vec[j]*w*params$theta
  }
  return(energy)
}

calculate_group_energy_A_c <- function(A_vec, L_vec, params_A, group_size=4) {
  energy <- 0
  for (i in 1:group_size) {
    int_nb <- get_interference_neighbors_c(i, group_size)
    energy <- energy + A_vec[i]*(params_A$gamma0 + params_A$gamma1*L_vec[i])
    if (length(int_nb)>0) energy <- energy + A_vec[i]*params_A$gamma2*mean(L_vec[int_nb])
  }
  for (pair in list(c(1,2), c(3,4))) {
    i <- pair[1]; j <- pair[2]
    w <- 0.5*(1/length(get_interference_neighbors_c(i,group_size)) +
                1/length(get_interference_neighbors_c(j,group_size)))
    energy <- energy + A_vec[i]*A_vec[j]*w*params_A$psi
  }
  return(energy)
}

calculate_group_energy_L_c <- function(L_vec, params_L, group_size=4) {
  energy <- sum(L_vec)*params_L$alpha
  for (i in 1:(group_size-1)) {
    j <- i+1
    dep_i <- get_dependence_neighbors_c(i, group_size)
    dep_j <- get_dependence_neighbors_c(j, group_size)
    w <- 0.5*(1/length(dep_i) + 1/length(dep_j))
    energy <- energy + L_vec[i]*L_vec[j]*w*params_L$omega
  }
  return(energy)
}

# -----------------------------------------------------------------------------
# Conditional probability functions — graph (c)
# -----------------------------------------------------------------------------

prob_Y_given_rest_c <- function(i, Y, A, L, params, network) {
  gs <- network$group_size; g <- ceiling(i/gs); il <- i-(g-1)*gs
  idx <- ((g-1)*gs+1):(g*gs)
  Yg <- Y[idx]; Ag <- A[idx]; Lg <- L[idx]
  Y1 <- Yg; Y1[il] <- 1; Y0 <- Yg; Y0[il] <- 0
  e1 <- calculate_group_energy_Y_c(Y1, Ag, Lg, params, gs)
  e0 <- calculate_group_energy_Y_c(Y0, Ag, Lg, params, gs)
  me <- max(e0,e1); exp(e1-me)/(exp(e0-me)+exp(e1-me))
}

prob_A_given_rest_c <- function(i, A, L, params_A, network) {
  gs <- network$group_size; g <- ceiling(i/gs); il <- i-(g-1)*gs
  idx <- ((g-1)*gs+1):(g*gs)
  Ag <- A[idx]; Lg <- L[idx]
  A1 <- Ag; A1[il] <- 1; A0 <- Ag; A0[il] <- 0
  e1 <- calculate_group_energy_A_c(A1, Lg, params_A, gs)
  e0 <- calculate_group_energy_A_c(A0, Lg, params_A, gs)
  me <- max(e0,e1); exp(e1-me)/(exp(e0-me)+exp(e1-me))
}

prob_L_given_rest_c <- function(i, L, params_L, network) {
  gs <- network$group_size; g <- ceiling(i/gs); il <- i-(g-1)*gs
  idx <- ((g-1)*gs+1):(g*gs)
  Lg <- L[idx]
  L1 <- Lg; L1[il] <- 1; L0 <- Lg; L0[il] <- 0
  e1 <- calculate_group_energy_L_c(L1, params_L, gs)
  e0 <- calculate_group_energy_L_c(L0, params_L, gs)
  me <- max(e0,e1); exp(e1-me)/(exp(e0-me)+exp(e1-me))
}

# -----------------------------------------------------------------------------
# Analytical estimands — graph (c)
# -----------------------------------------------------------------------------

calculate_psi_analytical_group_c <- function(i, a_vec, params, params_L, group_size=4) {
  L_cfg <- generate_binary_configs(group_size)
  Y_cfg <- generate_binary_configs(group_size)
  fL    <- exp(apply(L_cfg,1,function(l) calculate_group_energy_L_c(l,params_L,group_size)))
  fL    <- fL/sum(fL)
  EYi   <- sapply(1:nrow(L_cfg), function(li) {
    fY <- exp(apply(Y_cfg,1,function(y) calculate_group_energy_Y_c(y,a_vec,L_cfg[li,],params,group_size)))
    sum(Y_cfg[,i]*fY/sum(fY))
  })
  return(sum(EYi*fL))
}

calculate_prob_A_given_L_c <- function(A_vec, L_vec, params_A, group_size=4) {
  A_cfg <- generate_binary_configs(group_size)
  fA    <- exp(apply(A_cfg,1,function(a) calculate_group_energy_A_c(a,L_vec,params_A,group_size)))
  a_idx <- which(apply(A_cfg,1,function(x) all(x==A_vec)))
  return(fA[a_idx]/sum(fA))
}

calculate_marginal_pi_A_c <- function(A_vec, params_A, params_L, group_size=4) {
  L_cfg <- generate_binary_configs(group_size)
  fL    <- exp(apply(L_cfg,1,function(l) calculate_group_energy_L_c(l,params_L,group_size)))
  fL    <- fL/sum(fL)
  sum(sapply(1:nrow(L_cfg), function(li)
    calculate_prob_A_given_L_c(A_vec,L_cfg[li,],params_A,group_size)*fL[li]))
}

calculate_estimands_analytical_c <- function(network, params, params_A, params_L) {
  gs    <- network$group_size; n_groups <- network$n_groups
  A_cfg <- generate_binary_configs(gs)
  cat("  Computing allocation probabilities pi(a)...\n")
  pi_A  <- sapply(1:nrow(A_cfg), function(k)
    calculate_marginal_pi_A_c(A_cfg[k,],params_A,params_L,gs))
  cat("  Allocation probabilities sum to:", round(sum(pi_A),6), "\n")
  psi_1 <- psi_0 <- psi_0_0 <- numeric(gs)
  cat("  Computing analytical psi for each position...\n")
  for (i in 1:gs) {
    cat("    Position", i, "\n")
    cfg1 <- which(A_cfg[,i]==1); pi1 <- pi_A[cfg1]/sum(pi_A[cfg1])
    psi_1[i] <- sum(pi1*sapply(cfg1,function(k)
      calculate_psi_analytical_group_c(i,A_cfg[k,],params,params_L,gs)))
    cfg0 <- which(A_cfg[,i]==0); pi0 <- pi_A[cfg0]/sum(pi_A[cfg0])
    psi_0[i] <- sum(pi0*sapply(cfg0,function(k)
      calculate_psi_analytical_group_c(i,A_cfg[k,],params,params_L,gs)))
    psi_0_0[i] <- calculate_psi_analytical_group_c(i,rep(0,gs),params,params_L,gs)
  }
  DE  <- mean(psi_1-psi_0)
  IE  <- mean(psi_0-psi_0_0)
  ATE <- mean(psi_1-psi_0_0)   # ATE = DE + IE
  return(list(DE=DE, IE=IE, ATE=ATE,
              psi_1=rep(psi_1,n_groups), psi_0=rep(psi_0,n_groups),
              psi_0_0=rep(psi_0_0,n_groups),
              psi_1_group=psi_1, psi_0_group=psi_0, psi_0_0_group=psi_0_0, pi_A=pi_A))
}

# -----------------------------------------------------------------------------
# Gibbs samplers — graph (c)
# -----------------------------------------------------------------------------

gibbs_sampler_observational_c <- function(network, params, params_A, params_L,
                                          n_iter=5000, burn_in=1000) {
  N <- network$N
  L <- rbinom(N,1,0.5); A <- rbinom(N,1,0.5); Y <- rbinom(N,1,0.5)
  n_s <- n_iter-burn_in
  Ls <- matrix(0,n_s,N); As <- matrix(0,n_s,N); Ys <- matrix(0,n_s,N); si <- 1
  for (m in 1:n_iter) {
    for (i in sample(1:N)) {
      L[i] <- rbinom(1,1,prob_L_given_rest_c(i,L,params_L,network))
      A[i] <- rbinom(1,1,prob_A_given_rest_c(i,A,L,params_A,network))
      Y[i] <- rbinom(1,1,prob_Y_given_rest_c(i,Y,A,L,params,network))
    }
    if (m>burn_in) { Ls[si,]<-L; As[si,]<-A; Ys[si,]<-Y; si<-si+1 }
  }
  return(list(L_samples=Ls, A_samples=As, Y_samples=Ys))
}

gibbs_sampler_interventional_c <- function(a_vec, network, params, params_L,
                                           n_iter=5000, burn_in=1000) {
  N <- network$N
  L <- rbinom(N,1,0.5); Y <- rbinom(N,1,0.5)
  n_s <- n_iter-burn_in
  Ls <- matrix(0,n_s,N); Ys <- matrix(0,n_s,N); si <- 1
  for (m in 1:n_iter) {
    for (i in sample(1:N)) {
      L[i] <- rbinom(1,1,prob_L_given_rest_c(i,L,params_L,network))
      Y[i] <- rbinom(1,1,prob_Y_given_rest_c(i,Y,a_vec,L,params,network))
    }
    if (m>burn_in) { Ls[si,]<-L; Ys[si,]<-Y; si<-si+1 }
  }
  return(list(L_samples=Ls, Y_samples=Ys))
}

# -----------------------------------------------------------------------------
# Single simulation run — graph (c), returns DE, IE, ATE
# -----------------------------------------------------------------------------

run_single_simulation_c <- function(sim_id, network, params, params_A, params_L,
                                    n_iter, burn_in, n_alloc, seed=NULL) {
  if (!is.null(seed)) set.seed(seed+sim_id)
  N   <- network$N
  obs <- gibbs_sampler_observational_c(network, params, params_A, params_L,
                                       n_iter=n_alloc*50, burn_in=500)
  aidx  <- seq(1, nrow(obs$A_samples), length.out=n_alloc)
  psi_1 <- psi_0 <- cnt1 <- cnt0 <- numeric(N)
  for (k in 1:n_alloc) {
    a_vec  <- obs$A_samples[round(aidx[k]),]
    Ymeans <- colMeans(gibbs_sampler_interventional_c(a_vec,network,params,params_L,
                                                      n_iter,burn_in)$Y_samples)
    for (i in 1:N) {
      if (a_vec[i]==1) { psi_1[i]<-psi_1[i]+Ymeans[i]; cnt1[i]<-cnt1[i]+1 }
      else             { psi_0[i]<-psi_0[i]+Ymeans[i]; cnt0[i]<-cnt0[i]+1 }
    }
  }
  psi_1   <- psi_1/pmax(cnt1,1)
  psi_0   <- psi_0/pmax(cnt0,1)
  psi_0_0 <- colMeans(gibbs_sampler_interventional_c(rep(0,N),network,params,params_L,
                                                     n_iter,burn_in)$Y_samples)
  DE  <- mean(psi_1-psi_0)
  IE  <- mean(psi_0-psi_0_0)
  ATE <- mean(psi_1-psi_0_0)   # ATE = DE + IE
  return(list(DE=DE, IE=IE, ATE=ATE))
}

# -----------------------------------------------------------------------------
# Main parallel simulation — graph (c)
# -----------------------------------------------------------------------------

run_simulation_c_parallel <- function(
    n_sim=500, N=20,
    params   = list(beta0=-1, beta1=0.5, beta2=0.3, beta3=0.4, beta4=0.2, theta=0.3),
    params_A = list(gamma0=0, gamma1=0.3, gamma2=0.2, psi=0.4),
    params_L = list(alpha=-0.5, omega=0.4),
    n_iter=5000, burn_in=1000, n_alloc=50,
    n_cores=NULL, verbose=TRUE, seed=42) {
  
  if (is.null(n_cores)) n_cores <- max(1, detectCores()-1)
  network <- create_network_structure_c(N)
  
  if (verbose) {
    cat("===========================================\n")
    cat("PARALLEL Simulation for Chain Graph (c)\n")
    cat("Structure: Interference = pairs (1-2),(3-4)\n")
    cat("           Dependence   = full chain 1-2-3-4\n")
    cat("N =", N, "| n_sim =", n_sim, "| n_cores =", n_cores, "\n")
    cat("===========================================\n\n")
  }
  
  if (verbose) cat("Calculating analytical truth...\n")
  truth <- calculate_estimands_analytical_c(network, params, params_A, params_L)
  if (verbose) cat("\nTrue DE:", round(truth$DE,4),
                   "| True IE:", round(truth$IE,4),
                   "| True ATE:", round(truth$ATE,4), "\n\n")
  
  cl <- makeCluster(n_cores); registerDoParallel(cl)
  clusterExport(cl, c(
    "create_network_structure_c", "generate_binary_configs",
    "get_interference_neighbors_c", "get_dependence_neighbors_c",
    "calculate_group_energy_Y_c", "calculate_group_energy_A_c", "calculate_group_energy_L_c",
    "prob_Y_given_rest_c", "prob_A_given_rest_c", "prob_L_given_rest_c",
    "gibbs_sampler_observational_c", "gibbs_sampler_interventional_c",
    "run_single_simulation_c"), envir=environment())
  
  if (verbose) cat("Running", n_sim, "simulations in parallel...\n")
  t0 <- Sys.time()
  res_list <- foreach(sim=1:n_sim, .combine=rbind) %dopar% {
    r <- run_single_simulation_c(sim, network, params, params_A, params_L,
                                 n_iter, burn_in, n_alloc, seed)
    c(DE=r$DE, IE=r$IE, ATE=r$ATE)
  }
  stopCluster(cl)
  elapsed <- difftime(Sys.time(), t0, units="mins")
  if (verbose) cat("Completed in", round(as.numeric(elapsed),2), "minutes\n\n")
  
  DE_est  <- res_list[,"DE"]
  IE_est  <- res_list[,"IE"]
  ATE_est <- res_list[,"ATE"]
  
  coverage <- function(est, tv) mean(abs(est-tv) <= 1.96*sd(est))
  
  results <- list(
    true_DE=truth$DE, true_IE=truth$IE, true_ATE=truth$ATE,
    DE_estimates=DE_est, IE_estimates=IE_est, ATE_estimates=ATE_est,
    DE_mean=mean(DE_est),   DE_bias=mean(DE_est)-truth$DE,
    DE_se=sd(DE_est),       DE_rmse=sqrt(mean((DE_est-truth$DE)^2)),
    IE_mean=mean(IE_est),   IE_bias=mean(IE_est)-truth$IE,
    IE_se=sd(IE_est),       IE_rmse=sqrt(mean((IE_est-truth$IE)^2)),
    ATE_mean=mean(ATE_est), ATE_bias=mean(ATE_est)-truth$ATE,
    ATE_se=sd(ATE_est),     ATE_rmse=sqrt(mean((ATE_est-truth$ATE)^2)),
    DE_coverage=coverage(DE_est,truth$DE),
    IE_coverage=coverage(IE_est,truth$IE),
    ATE_coverage=coverage(ATE_est,truth$ATE),
    elapsed_time_mins=as.numeric(elapsed), n_cores_used=n_cores)
  
  if (verbose) {
    cat("===========================================\n")
    cat("RESULTS FOR GRAPH (c)\n")
    cat("===========================================\n")
    for (nm in c("DE","IE","ATE"))
      cat(sprintf("\n%s:\n  True: %.4f | Mean: %.4f | Bias: %.4f | SE: %.4f | RMSE: %.4f | 95%% Cov: %.1f%%\n",
                  nm, results[[paste0("true_",nm)]], results[[paste0(nm,"_mean")]],
                  results[[paste0(nm,"_bias")]], results[[paste0(nm,"_se")]],
                  results[[paste0(nm,"_rmse")]], results[[paste0(nm,"_coverage")]]*100))
    cat("\nComputation time:", round(results$elapsed_time_mins,2), "minutes\n")
    cat("Cores used:", results$n_cores_used, "\n")
  }
  return(results)
}


# =============================================================================
# RUN ALL THREE GRAPHS
# =============================================================================

params   <- list(beta0=-1, beta1=0.5, beta2=0.3, beta3=0.4, beta4=0.2, theta=0.3)
params_A <- list(gamma0=0, gamma1=0.3, gamma2=0.2, psi=0.4)
params_L <- list(alpha=-0.5, omega=0.4)

cat("Detected", detectCores(), "CPU cores\n\n")

results_a <- run_simulation_a_parallel(
  n_sim=500, N=20, params=params, params_A=params_A, params_L=params_L,
  n_iter=5000, burn_in=1000, n_alloc=50, n_cores=NULL, verbose=TRUE, seed=42)

results_b <- run_simulation_b_parallel(
  n_sim=500, N=20, params=params, params_A=params_A, params_L=params_L,
  n_iter=5000, burn_in=1000, n_alloc=50, n_cores=NULL, verbose=TRUE, seed=42)

results_c <- run_simulation_c_parallel(
  n_sim=500, N=20, params=params, params_A=params_A, params_L=params_L,
  n_iter=5000, burn_in=1000, n_alloc=50, n_cores=NULL, verbose=TRUE, seed=42)

