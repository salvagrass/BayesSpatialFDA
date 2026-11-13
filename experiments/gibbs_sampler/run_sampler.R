# experiments/gibbs_sampler/run_sampler.R
# Main wrapper for the Marginalised Adaptive Anisotropic SPDE Gibbs Sampler

#' Run the Marginalised Adaptive Gibbs Sampler (R Prototype)
#'
#' @param Y_mat Numeric matrix (N x M) of observations across N vertices and M years
#' @param matrices List containing SPDE FEM matrices (`C`, `K11`, `K22`, `K12`)
#' @param priors List of prior hyperparameters (`ig_a`, `ig_b`, `mu_mean`, `mu_prec`, `hyper_sd`)
#' @param control List of MCMC settings (`G`, `burnin`, `thin`, `adapt_start`, `adapt_every`, `init_sd`, `scale_shrink`, `progress`)
#' @param init Optional list of initial values (`mu`, `beta`, `sigma2`, `tau2`, `kappa2`, `gamma`, `theta`)
#' @return A structured list containing posterior draws, online summaries, and diagnostics
run_gibbs_sampler <- function(
  Y_mat,
  matrices,
  priors = list(),
  control = list(),
  init = list()
) {
  N <- nrow(Y_mat)
  M <- ncol(Y_mat)

  # --- 1. Unpack Control Settings ---
  G <- control$G %||% 5000L
  burnin <- control$burnin %||% 1000L
  thin <- control$thin %||% 4L
  adapt_start <- control$adapt_start %||% 200L
  adapt_every <- control$adapt_every %||% 10L
  init_sd <- control$init_sd %||% 0.05
  scale_shrink <- control$scale_shrink %||% 0.7
  show_prog <- control$progress %||% TRUE

  if (burnin >= G) {
    stop("'burnin' must be strictly less than 'G'.")
  }
  idx_post <- (burnin + 1L):G
  n_post <- length(idx_post)

  # --- 2. Unpack Priors ---
  prior_ig_a <- priors$ig_a %||% 0.01
  prior_ig_b <- priors$ig_b %||% 0.01
  prior_mu_mean <- priors$mu_mean %||% rep(0, N)
  prior_mu_prec <- priors$mu_prec %||% Matrix::Diagonal(N, 1 / 100)
  hyper_sd <- priors$hyper_sd %||% 5.0

  # --- 3. Extract FEM Matrices ---
  C_mat <- matrices$C
  K11 <- matrices$K11
  K22 <- matrices$K22
  K12 <- matrices$K12
  inv_diag_C <- 1 / Matrix::diag(C_mat)

  # --- 4. Initialise States ---
  curr_mu <- init$mu %||% rep(0, N)
  curr_beta <- init$beta %||% matrix(0, nrow = N, ncol = M)
  curr_sigma2 <- init$sigma2 %||% 1.0
  curr_tau2 <- init$tau2 %||% 0.5
  curr_kappa2 <- init$kappa2 %||% 0.5
  curr_gamma <- init$gamma %||% 1.0
  curr_theta <- init$theta %||% 0.0

  Q_pr_curr <- build_spde_precision(
    tau2 = curr_tau2,
    kappa2 = curr_kappa2,
    gamma = curr_gamma,
    theta = curr_theta,
    C_mat = C_mat,
    K11 = K11,
    K22 = K22,
    K12 = K12,
    inv_diag_C = inv_diag_C
  )
  chol_pr_curr <- Matrix::chol(Q_pr_curr)

  # --- 5. Storage & Adaptive Proposal Setup ---
  samples_sigma2 <- numeric(G)
  samples_tau2 <- numeric(G)
  samples_kappa2 <- numeric(G)
  samples_gamma <- numeric(G)
  samples_theta <- numeric(G)
  samples_beta1 <- numeric(G)

  # Online sufficient statistics for beta and mu (post burn-in)
  beta_sum_post <- matrix(0, nrow = N, ncol = M)
  beta_sq_sum_post <- matrix(0, nrow = N, ncol = M)
  mu_sum_post <- numeric(N)
  mu_sq_sum_post <- numeric(N)

  # Thinned storage for post burn-in beta draws
  keep_iter <- idx_post[idx_post %% thin == 0L]
  beta_draws <- array(NA_real_, dim = c(N, M, length(keep_iter)))
  keep_idx <- 1L

  # Adaptive random-walk covariance
  mcmc_history <- matrix(0, nrow = G, ncol = 4)
  colnames(mcmc_history) <- c("log_tau2", "log_kappa2", "log_gamma", "theta")
  prop_cov <- diag(x = init_sd^2, nrow = 4)
  scale_factor <- ((2.38^2) / 4) * scale_shrink
  acc_count <- 0L

  # Precompute constant posterior shape for sigma2
  a_post_sig <- (N * M) / 2 + prior_ig_a

  if (show_prog) {
    cat(sprintf(
      "Starting Marginalised Adaptive Gibbs Sampler (N = %d, M = %d, G = %d)...\n",
      N,
      M,
      G
    ))
    pb <- txtProgressBar(min = 0, max = G, style = 3)
  }

  start_time <- Sys.time()

  # --- 6. Main MCMC Loop ---
  for (g in seq_len(G)) {
    # Step 1: Gibbs update of long-term mean mu
    curr_mu <- update_mu(
      N = N,
      M = M,
      prior_mean = prior_mu_mean,
      prior_prec = prior_mu_prec,
      Y_mat = Y_mat,
      curr_beta = curr_beta,
      curr_sigma2 = curr_sigma2
    )

    # Step 2: Gibbs update of observation noise variance sigma2
    eps_mat <- Y_mat - curr_mu - curr_beta
    b_post_sig <- sum(eps_mat^2) / 2 + prior_ig_b
    curr_sigma2 <- 1 / rgamma(1, shape = a_post_sig, rate = b_post_sig)

    # Step 3: Marginalised M-H update of SPDE hyperparameters
    Y_star <- Y_mat - curr_mu
    B_mat <- Y_star / curr_sigma2
    Q_post_curr <- Q_pr_curr + Matrix::Diagonal(N, 1 / curr_sigma2)
    chol_post_curr <- Matrix::chol(Q_post_curr)
    M_post_curr <- Matrix::solve(
      chol_post_curr,
      Matrix::solve(t(chol_post_curr), B_mat)
    )

    log_det_pr_curr <- 2 * sum(log(Matrix::diag(chol_pr_curr)))
    log_det_post_curr <- 2 * sum(log(Matrix::diag(chol_post_curr)))
    quad_curr <- sum(B_mat * M_post_curr)
    log_lik_curr <- 0.5 *
      M *
      log_det_pr_curr -
      0.5 * M * log_det_post_curr +
      0.5 * quad_curr

    if (g > adapt_start && (g %% adapt_every == 0L)) {
      emp_cov <- cov(mcmc_history[1:(g - 1L), , drop = FALSE])
      prop_cov <- scale_factor * emp_cov + diag(1e-4, nrow = 4)
    }

    mh_step <- marginal_metropolis_hastings(
      curr_tau2 = curr_tau2,
      curr_kappa2 = curr_kappa2,
      curr_gamma = curr_gamma,
      curr_theta = curr_theta,
      B_mat = B_mat,
      curr_sigma2 = curr_sigma2,
      C_mat = C_mat,
      K11 = K11,
      K22 = K22,
      K12 = K12,
      inv_diag_C = inv_diag_C,
      prop_cov = prop_cov,
      Q_pr_curr = Q_pr_curr,
      chol_pr_curr = chol_pr_curr,
      chol_post_curr = chol_post_curr,
      M_post_curr = M_post_curr,
      log_lik_curr = log_lik_curr,
      prior_sd = hyper_sd
    )

    curr_tau2 <- mh_step$tau2
    curr_kappa2 <- mh_step$kappa2
    curr_gamma <- mh_step$gamma
    curr_theta <- mh_step$theta
    Q_pr_curr <- mh_step$Q_pr
    chol_pr_curr <- mh_step$chol_pr
    chol_post_curr <- mh_step$chol_post
    M_post_curr <- mh_step$M_post

    if (mh_step$accepted) {
      acc_count <- acc_count + 1L
    }

    # Step 4: Joint Gibbs draw of spatial anomaly field beta
    Z_mat <- matrix(rnorm(N * M), nrow = N, ncol = M)
    curr_beta <- as.matrix(M_post_curr + Matrix::solve(chol_post_curr, Z_mat))

    # --- Save Samples & Online Summaries ---
    samples_sigma2[g] <- curr_sigma2
    samples_tau2[g] <- curr_tau2
    samples_kappa2[g] <- curr_kappa2
    samples_gamma[g] <- curr_gamma
    samples_theta[g] <- curr_theta
    samples_beta1[g] <- curr_beta[1L, 1L]

    if (g > burnin) {
      beta_sum_post <- beta_sum_post + curr_beta
      beta_sq_sum_post <- beta_sq_sum_post + curr_beta^2
      mu_sum_post <- mu_sum_post + curr_mu
      mu_sq_sum_post <- mu_sq_sum_post + curr_mu^2

      if (keep_idx <= length(keep_iter) && g == keep_iter[keep_idx]) {
        beta_draws[,, keep_idx] <- curr_beta
        keep_idx <- keep_idx + 1L
      }
    }

    mcmc_history[g, ] <- c(
      log(curr_tau2),
      log(curr_kappa2),
      log(curr_gamma),
      curr_theta
    )

    if (show_prog) setTxtProgressBar(pb, g)
  }

  if (show_prog) {
    close(pb)
  }
  elapsed_sec <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # --- 7. Compute Post Burn-In Moments ---
  beta_mean <- beta_sum_post / n_post
  beta_sd <- sqrt(pmax((beta_sq_sum_post / n_post) - beta_mean^2, 0))
  mu_mean <- mu_sum_post / n_post
  mu_sd <- sqrt(pmax(0, (mu_sq_sum_post / n_post) - mu_mean^2))

  return(list(
    draws = list(
      sigma2 = samples_sigma2,
      tau2 = samples_tau2,
      kappa2 = samples_kappa2,
      gamma = samples_gamma,
      theta = samples_theta,
      beta1 = samples_beta1,
      beta = beta_draws,
      keep_iter = keep_iter
    ),
    summary = list(
      beta_mean = beta_mean,
      beta_sd = beta_sd,
      mu_mean = mu_mean,
      mu_sd = mu_sd
    ),
    diagnostics = list(
      acc_rate = acc_count / G,
      elapsed_sec = elapsed_sec,
      G = G,
      burnin = burnin,
      thin = thin,
      idx_post = idx_post
    )
  ))
}
