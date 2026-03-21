# =============================================================================
# Revised Simulation Setting 2: Network Misspecification Study
# TRUE DGP: Graph (c)
# 20 units = 5 independent groups of size 4
# =============================================================================

set.seed(42)

# =============================================================================
# SECTION 1: GROUP-LEVEL NETWORK STRUCTURES (size = 4)
# =============================================================================

create_group_network_a <- function() {
  G_interference <- matrix(0, 4, 4)
  G_dependence   <- matrix(0, 4, 4)
  
  # chain interference: 1-2-3-4
  for (k in 1:3) {
    G_interference[k, k + 1] <- 1
    G_interference[k + 1, k] <- 1
  }
  
  # chain dependence: 1-2-3-4
  for (k in 1:3) {
    G_dependence[k, k + 1] <- 1
    G_dependence[k + 1, k] <- 1
  }
  
  list(G_interference = G_interference,
       G_dependence   = G_dependence,
       group_size     = 4)
}

create_group_network_b <- function() {
  G_interference <- matrix(0, 4, 4)
  G_dependence   <- matrix(0, 4, 4)
  
  # pair interference: (1,2), (3,4)
  G_interference[1, 2] <- G_interference[2, 1] <- 1
  G_interference[3, 4] <- G_interference[4, 3] <- 1
  
  # pair dependence: (1,2), (3,4)
  G_dependence[1, 2] <- G_dependence[2, 1] <- 1
  G_dependence[3, 4] <- G_dependence[4, 3] <- 1
  
  list(G_interference = G_interference,
       G_dependence   = G_dependence,
       group_size     = 4)
}

create_group_network_c <- function() {
  G_interference <- matrix(0, 4, 4)
  G_dependence   <- matrix(0, 4, 4)
  
  # pair interference: (1,2), (3,4)
  G_interference[1, 2] <- G_interference[2, 1] <- 1
  G_interference[3, 4] <- G_interference[4, 3] <- 1
  
  # chain dependence: 1-2-3-4
  for (k in 1:3) {
    G_dependence[k, k + 1] <- 1
    G_dependence[k + 1, k] <- 1
  }
  
  list(G_interference = G_interference,
       G_dependence   = G_dependence,
       group_size     = 4)
}

# =============================================================================
# SECTION 2: HELPER FUNCTIONS
# =============================================================================

generate_binary_configs <- function(n) {
  configs <- matrix(0, 2^n, n)
  for (j in 1:n) {
    configs[, j] <- rep(c(0, 1), each = 2^(n - j), times = 2^(j - 1))
  }
  configs
}

# =============================================================================
# SECTION 3: ENERGY FUNCTIONS
# =============================================================================

calculate_energy_L <- function(L_vec, params_L, network) {
  n <- length(L_vec)
  energy <- 0
  
  for (i in 1:n) {
    energy <- energy + L_vec[i] * params_L$alpha
  }
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (network$G_dependence[i, j] == 1) {
        # IMPORTANT: use unweighted edge contribution so missing edge 2-3 matters more
        energy <- energy + params_L$omega * L_vec[i] * L_vec[j]
      }
    }
  }
  
  energy
}

calculate_energy_A <- function(A_vec, L_vec, params_A, network) {
  n <- length(A_vec)
  energy <- 0
  
  for (i in 1:n) {
    neigh_i <- which(network$G_interference[i, ] == 1)
    
    linear_part <- params_A$gamma0 + params_A$gamma1 * L_vec[i]
    
    if (length(neigh_i) > 0) {
      linear_part <- linear_part + params_A$gamma2 * sum(L_vec[neigh_i])
    }
    
    energy <- energy + A_vec[i] * linear_part
  }
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (network$G_dependence[i, j] == 1) {
        energy <- energy + params_A$psi * A_vec[i] * A_vec[j]
      }
    }
  }
  
  energy
}

calculate_energy_Y <- function(Y_vec, A_vec, L_vec, params, network) {
  n <- length(Y_vec)
  energy <- 0
  
  for (i in 1:n) {
    neigh_i_interf <- which(network$G_interference[i, ] == 1)
    
    linear_part <- params$beta0 +
      params$beta1 * A_vec[i] +
      params$beta2 * L_vec[i]
    
    if (length(neigh_i_interf) > 0) {
      linear_part <- linear_part +
        params$beta3 * sum(A_vec[neigh_i_interf]) +
        params$beta4 * sum(L_vec[neigh_i_interf])
    }
    
    energy <- energy + Y_vec[i] * linear_part
  }
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (network$G_dependence[i, j] == 1) {
        energy <- energy + params$theta * Y_vec[i] * Y_vec[j]
      }
    }
  }
  
  energy
}

# =============================================================================
# SECTION 4: EXACT GROUP-LEVEL TRUTH UNDER TRUE MODEL (c)
# =============================================================================

calculate_f_L <- function(params_L, network_true) {
  L_configs <- generate_binary_configs(4)
  f_L_unnorm <- apply(L_configs, 1, function(l) {
    exp(calculate_energy_L(l, params_L, network_true))
  })
  f_L <- f_L_unnorm / sum(f_L_unnorm)
  list(L_configs = L_configs, probs = f_L)
}

calculate_f_A_given_L <- function(L_vec, params_A, network_true) {
  A_configs <- generate_binary_configs(4)
  f_A_unnorm <- apply(A_configs, 1, function(a) {
    exp(calculate_energy_A(a, L_vec, params_A, network_true))
  })
  f_A <- f_A_unnorm / sum(f_A_unnorm)
  list(A_configs = A_configs, probs = f_A)
}

calculate_EY_given_A_L <- function(A_vec, L_vec, params, network_true) {
  Y_configs <- generate_binary_configs(4)
  f_Y_unnorm <- apply(Y_configs, 1, function(y) {
    exp(calculate_energy_Y(y, A_vec, L_vec, params, network_true))
  })
  f_Y <- f_Y_unnorm / sum(f_Y_unnorm)
  colSums(Y_configs * f_Y)
}

# psi_i(a) = E[Yi(a)] for a specific intervention vector a
calculate_psi_vector_true <- function(A_vec, params, params_L, network_true) {
  fL_obj <- calculate_f_L(params_L, network_true)
  L_configs <- fL_obj$L_configs
  f_L <- fL_obj$probs
  
  psi <- rep(0, 4)
  for (l_idx in 1:nrow(L_configs)) {
    L_vec <- L_configs[l_idx, ]
    EY <- calculate_EY_given_A_L(A_vec, L_vec, params, network_true)
    psi <- psi + f_L[l_idx] * EY
  }
  psi
}

# marginal allocation distribution pi(a) = sum_l f(a|l) f(l)
calculate_pi_true <- function(params_A, params_L, network_true) {
  fL_obj <- calculate_f_L(params_L, network_true)
  L_configs <- fL_obj$L_configs
  f_L <- fL_obj$probs
  
  A_configs <- generate_binary_configs(4)
  pi_A <- rep(0, nrow(A_configs))
  
  for (l_idx in 1:nrow(L_configs)) {
    L_vec <- L_configs[l_idx, ]
    fA_obj <- calculate_f_A_given_L(L_vec, params_A, network_true)
    pi_A <- pi_A + f_L[l_idx] * fA_obj$probs
  }
  
  list(A_configs = A_configs, probs = pi_A)
}

calculate_true_estimands <- function(params, params_A, params_L) {
  network_true <- create_group_network_c()
  
  pi_obj <- calculate_pi_true(params_A, params_L, network_true)
  A_configs <- pi_obj$A_configs
  pi_A <- pi_obj$probs
  
  psi_1   <- rep(0, 4)
  psi_0   <- rep(0, 4)
  psi_00  <- calculate_psi_vector_true(rep(0, 4), params, params_L, network_true)
  
  for (i in 1:4) {
    idx1 <- which(A_configs[, i] == 1)
    idx0 <- which(A_configs[, i] == 0)
    
    w1 <- pi_A[idx1] / sum(pi_A[idx1])
    w0 <- pi_A[idx0] / sum(pi_A[idx0])
    
    tmp1 <- rep(0, 4)
    for (k in seq_along(idx1)) {
      tmp1 <- tmp1 + w1[k] * calculate_psi_vector_true(A_configs[idx1[k], ], params, params_L, network_true)
    }
    psi_1[i] <- tmp1[i]
    
    tmp0 <- rep(0, 4)
    for (k in seq_along(idx0)) {
      tmp0 <- tmp0 + w0[k] * calculate_psi_vector_true(A_configs[idx0[k], ], params, params_L, network_true)
    }
    psi_0[i] <- tmp0[i]
  }
  
  DE_unit <- psi_1 - psi_0
  IE_unit <- psi_0 - psi_00
  
  list(
    DE = mean(DE_unit),
    IE = mean(IE_unit),
    ATE = mean(psi_1 - psi_00),
    DE_unit = DE_unit,
    IE_unit = IE_unit,
    psi_1 = psi_1,
    psi_0 = psi_0,
    psi_00 = psi_00
  )
}

# =============================================================================
# SECTION 5: GIBBS SAMPLERS AT GROUP LEVEL
# =============================================================================

gibbs_sample_group_true_c <- function(params, params_A, params_L,
                                      n_iter = 4000, burn_in = 1000) {
  network_true <- create_group_network_c()
  n <- 4
  
  L <- rbinom(n, 1, 0.5)
  A <- rbinom(n, 1, 0.5)
  Y <- rbinom(n, 1, 0.5)
  
  n_keep <- n_iter - burn_in
  L_samples <- matrix(0, n_keep, n)
  A_samples <- matrix(0, n_keep, n)
  Y_samples <- matrix(0, n_keep, n)
  
  keep_idx <- 1
  
  for (m in 1:n_iter) {
    for (i in sample(1:n)) {
      # update L_i
      L1 <- L; L1[i] <- 1
      L0 <- L; L0[i] <- 0
      e1 <- calculate_energy_L(L1, params_L, network_true)
      e0 <- calculate_energy_L(L0, params_L, network_true)
      p1 <- exp(e1) / (exp(e1) + exp(e0))
      L[i] <- rbinom(1, 1, p1)
      
      # update A_i
      A1 <- A; A1[i] <- 1
      A0 <- A; A0[i] <- 0
      e1 <- calculate_energy_A(A1, L, params_A, network_true)
      e0 <- calculate_energy_A(A0, L, params_A, network_true)
      p1 <- exp(e1) / (exp(e1) + exp(e0))
      A[i] <- rbinom(1, 1, p1)
      
      # update Y_i
      Y1 <- Y; Y1[i] <- 1
      Y0 <- Y; Y0[i] <- 0
      e1 <- calculate_energy_Y(Y1, A, L, params, network_true)
      e0 <- calculate_energy_Y(Y0, A, L, params, network_true)
      p1 <- exp(e1) / (exp(e1) + exp(e0))
      Y[i] <- rbinom(1, 1, p1)
    }
    
    if (m > burn_in) {
      L_samples[keep_idx, ] <- L
      A_samples[keep_idx, ] <- A
      Y_samples[keep_idx, ] <- Y
      keep_idx <- keep_idx + 1
    }
  }
  
  list(L_samples = L_samples, A_samples = A_samples, Y_samples = Y_samples)
}

# Under fitted graph, estimate E[Y(a)] by integrating over L with Gibbs
gibbs_intervene_group <- function(a_vec, params, params_L, network_fit,
                                  n_iter = 4000, burn_in = 1000) {
  n <- 4
  L <- rbinom(n, 1, 0.5)
  Y <- rbinom(n, 1, 0.5)
  
  n_keep <- n_iter - burn_in
  Y_samples <- matrix(0, n_keep, n)
  
  keep_idx <- 1
  
  for (m in 1:n_iter) {
    for (i in sample(1:n)) {
      # update L_i
      L1 <- L; L1[i] <- 1
      L0 <- L; L0[i] <- 0
      e1 <- calculate_energy_L(L1, params_L, network_fit)
      e0 <- calculate_energy_L(L0, params_L, network_fit)
      p1 <- exp(e1) / (exp(e1) + exp(e0))
      L[i] <- rbinom(1, 1, p1)
      
      # update Y_i under intervention A = a_vec
      Y1 <- Y; Y1[i] <- 1
      Y0 <- Y; Y0[i] <- 0
      e1 <- calculate_energy_Y(Y1, a_vec, L, params, network_fit)
      e0 <- calculate_energy_Y(Y0, a_vec, L, params, network_fit)
      p1 <- exp(e1) / (exp(e1) + exp(e0))
      Y[i] <- rbinom(1, 1, p1)
    }
    
    if (m > burn_in) {
      Y_samples[keep_idx, ] <- Y
      keep_idx <- keep_idx + 1
    }
  }
  
  colMeans(Y_samples)
}

# =============================================================================
# SECTION 6: ONE 20-UNIT REPLICATE
# =============================================================================

estimate_one_replicate <- function(params, params_A, params_L,
                                   network_fit,
                                   n_groups = 5,
                                   n_alloc = 40,
                                   n_iter_data = 3000,
                                   burn_in_data = 1000,
                                   n_iter_intervene = 3000,
                                   burn_in_intervene = 1000) {
  
  # STEP 1: Generate observed A allocations from true model (c)
  # We do this group by group, then pool the units
  A_obs_all_groups <- array(0, dim = c(n_alloc, n_groups, 4))
  
  for (g in 1:n_groups) {
    obs <- gibbs_sample_group_true_c(params, params_A, params_L,
                                     n_iter = n_alloc * 60,
                                     burn_in = 1000)
    pick <- round(seq(1, nrow(obs$A_samples), length.out = n_alloc))
    A_obs_all_groups[, g, ] <- obs$A_samples[pick, ]
  }
  
  # STEP 2: For each unit, estimate psi_i(1), psi_i(0)
  psi_1_sum <- rep(0, n_groups * 4)
  psi_0_sum <- rep(0, n_groups * 4)
  count_1   <- rep(0, n_groups * 4)
  count_0   <- rep(0, n_groups * 4)
  
  for (k in 1:n_alloc) {
    for (g in 1:n_groups) {
      a_vec <- A_obs_all_groups[k, g, ]
      
      EY <- gibbs_intervene_group(a_vec, params, params_L, network_fit,
                                  n_iter = n_iter_intervene,
                                  burn_in = burn_in_intervene)
      
      idx_global <- ((g - 1) * 4 + 1):(g * 4)
      
      for (i in 1:4) {
        if (a_vec[i] == 1) {
          psi_1_sum[idx_global[i]] <- psi_1_sum[idx_global[i]] + EY[i]
          count_1[idx_global[i]]   <- count_1[idx_global[i]] + 1
        } else {
          psi_0_sum[idx_global[i]] <- psi_0_sum[idx_global[i]] + EY[i]
          count_0[idx_global[i]]   <- count_0[idx_global[i]] + 1
        }
      }
    }
  }
  
  psi_1_hat <- psi_1_sum / pmax(count_1, 1)
  psi_0_hat <- psi_0_sum / pmax(count_0, 1)
  
  # STEP 3: Estimate psi_i(0,0,0,0) in each group
  psi_00_hat <- rep(0, n_groups * 4)
  a_zero <- rep(0, 4)
  
  for (g in 1:n_groups) {
    EY0 <- gibbs_intervene_group(a_zero, params, params_L, network_fit,
                                 n_iter = n_iter_intervene,
                                 burn_in = burn_in_intervene)
    idx_global <- ((g - 1) * 4 + 1):(g * 4)
    psi_00_hat[idx_global] <- EY0
  }
  
  DE_hat  <- mean(psi_1_hat - psi_0_hat)
  IE_hat  <- mean(psi_0_hat - psi_00_hat)
  ATE_hat <- mean(psi_1_hat - psi_00_hat)
  
  # also report specifically for units 2 and 3 in each group
  idx23 <- c()
  for (g in 1:n_groups) {
    idx23 <- c(idx23, (g - 1) * 4 + 2, (g - 1) * 4 + 3)
  }
  
  DE_23_hat  <- mean((psi_1_hat - psi_0_hat)[idx23])
  IE_23_hat  <- mean((psi_0_hat - psi_00_hat)[idx23])
  ATE_23_hat <- mean((psi_1_hat - psi_00_hat)[idx23])
  
  list(
    DE = DE_hat,
    IE = IE_hat,
    ATE = ATE_hat,
    DE_23 = DE_23_hat,
    IE_23 = IE_23_hat,
    ATE_23 = ATE_23_hat
  )
}

# =============================================================================
# SECTION 7: MAIN SIMULATION
# =============================================================================

summarize_metric <- function(estimates, truth) {
  c(
    mean = mean(estimates),
    bias = mean(estimates) - truth,
    emp_sd = sd(estimates),
    rmse = sqrt(mean((estimates - truth)^2))
  )
}

run_simulation_setting2_revised <- function(
    n_sim = 200,
    n_groups = 5,
    params = list(beta0 = -1, beta1 = 0.5, beta2 = 0.3,
                  beta3 = 0.7, beta4 = 0.4, theta = 1.0),
    params_A = list(gamma0 = 0, gamma1 = 0.3, gamma2 = 0.3, psi = 0.6),
    params_L = list(alpha = -0.5, omega = 1.0),
    n_alloc = 40,
    n_iter_intervene = 2500,
    burn_in_intervene = 800,
    verbose = TRUE
) {
  
  net_a <- create_group_network_a()
  net_b <- create_group_network_b()
  net_c <- create_group_network_c()
  
  truth_group <- calculate_true_estimands(params, params_A, params_L)
  
  if (verbose) {
    cat("==============================================\n")
    cat("REVISED SETTING 2: NETWORK MISSPECIFICATION\n")
    cat("==============================================\n")
    cat("Total units:", n_groups * 4, "\n")
    cat("Number of groups:", n_groups, "\n")
    cat("True group-level graph: (c)\n")
    cat("  Interference: (1,2), (3,4)\n")
    cat("  Dependence:   (1,2), (2,3), (3,4)\n\n")
    cat("True DE :", round(truth_group$DE, 4), "\n")
    cat("True IE :", round(truth_group$IE, 4), "\n")
    cat("True ATE:", round(truth_group$ATE, 4), "\n\n")
    cat("Unit-specific truth within a group:\n")
    cat("DE by unit :", paste(round(truth_group$DE_unit, 4), collapse = ", "), "\n")
    cat("IE by unit :", paste(round(truth_group$IE_unit, 4), collapse = ", "), "\n")
    cat("==============================================\n\n")
  }
  
  # storage
  res <- list(
    a = matrix(NA, n_sim, 6),
    b = matrix(NA, n_sim, 6),
    c = matrix(NA, n_sim, 6)
  )
  colnames(res$a) <- colnames(res$b) <- colnames(res$c) <-
    c("DE", "IE", "ATE", "DE_23", "IE_23", "ATE_23")
  
  for (s in 1:n_sim) {
    if (verbose && s %% 20 == 0) {
      cat("Simulation", s, "of", n_sim, "\n")
    }
    
    est_a <- estimate_one_replicate(
      params, params_A, params_L, network_fit = net_a,
      n_groups = n_groups, n_alloc = n_alloc,
      n_iter_intervene = n_iter_intervene,
      burn_in_intervene = burn_in_intervene
    )
    
    est_b <- estimate_one_replicate(
      params, params_A, params_L, network_fit = net_b,
      n_groups = n_groups, n_alloc = n_alloc,
      n_iter_intervene = n_iter_intervene,
      burn_in_intervene = burn_in_intervene
    )
    
    est_c <- estimate_one_replicate(
      params, params_A, params_L, network_fit = net_c,
      n_groups = n_groups, n_alloc = n_alloc,
      n_iter_intervene = n_iter_intervene,
      burn_in_intervene = burn_in_intervene
    )
    
    res$a[s, ] <- unlist(est_a)
    res$b[s, ] <- unlist(est_b)
    res$c[s, ] <- unlist(est_c)
  }
  
  summary_table <- rbind(
    cbind(model = "Graph (a)", estimand = "DE",
          t(summarize_metric(res$a[, "DE"], truth_group$DE))),
    cbind(model = "Graph (a)", estimand = "IE",
          t(summarize_metric(res$a[, "IE"], truth_group$IE))),
    cbind(model = "Graph (a)", estimand = "ATE",
          t(summarize_metric(res$a[, "ATE"], truth_group$ATE))),
    
    cbind(model = "Graph (b)", estimand = "DE",
          t(summarize_metric(res$b[, "DE"], truth_group$DE))),
    cbind(model = "Graph (b)", estimand = "IE",
          t(summarize_metric(res$b[, "IE"], truth_group$IE))),
    cbind(model = "Graph (b)", estimand = "ATE",
          t(summarize_metric(res$b[, "ATE"], truth_group$ATE))),
    
    cbind(model = "Graph (c)", estimand = "DE",
          t(summarize_metric(res$c[, "DE"], truth_group$DE))),
    cbind(model = "Graph (c)", estimand = "IE",
          t(summarize_metric(res$c[, "IE"], truth_group$IE))),
    cbind(model = "Graph (c)", estimand = "ATE",
          t(summarize_metric(res$c[, "ATE"], truth_group$ATE)))
  )
  
  summary_table <- as.data.frame(summary_table, stringsAsFactors = FALSE)
  
  for (nm in c("mean", "bias", "emp_sd", "rmse")) {
    summary_table[[nm]] <- as.numeric(summary_table[[nm]])
  }
  
  list(
    truth = truth_group,
    raw = res,
    summary = summary_table
  )
}

# =============================================================================
# SECTION 8: RUN
# =============================================================================

params <- list(
  beta0 = -1,
  beta1 = 0.5,
  beta2 = 0.3,
  beta3 = 0.7,   # stronger interference effect
  beta4 = 0.4,
  theta = 1.0    # stronger Y-dependence
)

params_A <- list(
  gamma0 = 0,
  gamma1 = 0.3,
  gamma2 = 0.3,
  psi = 0.6
)

params_L <- list(
  alpha = -0.5,
  omega = 1.0    # stronger L-dependence
)

results_revised <- run_simulation_setting2_revised(
  n_sim = 100,              # increase to 300-500 for final table
  n_groups = 5,
  params = params,
  params_A = params_A,
  params_L = params_L,
  n_alloc = 40,
  n_iter_intervene = 2000,
  burn_in_intervene = 600,
  verbose = TRUE
)

print(results_revised$summary)
