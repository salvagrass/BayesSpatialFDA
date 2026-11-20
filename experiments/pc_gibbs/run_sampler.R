run_gibbs_sampler <- function(
  N,
  V,
  M,
  Phi,
  Y_mat,
  init_state,
  priors,
  matrices,
  num_iterations,
  num_burn_in,
  standard_step,
  verbose = FALSE
) {
  stopifnot(
    num_iterations > num_burn_in,
    nrow(Phi) == N,
    ncol(Phi) == V,
    nrow(Y_mat) == N,
    ncol(Y_mat) == M,
    length(init_state$mu) == V,
    identical(dim(init_state$beta_mat), c(V, M))
  )

  num_saved <- num_iterations - num_burn_in

  # Setting up the samples (latent fields live on the V mesh vertices)
  beta_samples <- array(data = NA_real_, dim = c(V, M, num_saved))
  mu_samples <- matrix(data = NA_real_, nrow = V, ncol = num_saved)
  theta_samples <- matrix(
    data = NA_real_,
    nrow = 4,
    ncol = num_iterations,
    dimnames = list(c("kappa", "sigma_beta", "v1", "v2"), NULL)
  )
  sigma2_samples <- numeric(num_iterations)

  num_accepted <- 0L
  curr_state <- init_state

  # Setup Adaptive Metropolis-Hastings (Logarithmic space)
  num_adaptive <- floor(0.3 * num_burn_in)
  eta_mean <- numeric(4)
  eta_cov_mat <- matrix(0, nrow = 4, ncol = 4)
  opt_rate <- (2.38^2) / 4
  step_rate <- standard_step

  if (verbose) {
    pb <- txtProgressBar(min = 0, max = num_iterations, style = 3)
    on.exit(close(pb), add = TRUE)
  }

  for (i in seq_len(num_iterations)) {
    # 1. Updating sigma2_epsilon
    curr_state$sigma2_epsilon <- update_sigma2_epsilon(
      N,
      V,
      M,
      Phi,
      Y_mat,
      curr_state$beta_mat,
      priors
    )

    # 2. First step for mu (centered)
    mu_half <- update_mu_centered(
      V,
      M,
      curr_state$beta_mat,
      curr_state$theta,
      matrices,
      priors
    )

    # 3. ASIS step (non-centered) updating both beta and mu
    res_asis <- asis_step(
      V,
      M,
      Phi,
      Y_mat,
      curr_state$beta_mat,
      mu_half,
      curr_state$sigma2_epsilon,
      priors
    )
    curr_state$mu <- res_asis$mu
    curr_state$beta_mat <- res_asis$beta_mat

    # 4. Metropolis-Hastings step for theta + collapsed beta update
    res_mh <- metropolis_step(
      N,
      V,
      M,
      Phi,
      step_rate,
      curr_state$theta,
      curr_state$mu,
      Y_mat,
      curr_state$sigma2_epsilon,
      priors,
      matrices
    )
    curr_state$theta <- res_mh$spde_theta
    curr_state$beta_mat <- res_mh$beta_mat

    # 5. Adattamento della covarianza della proposta durante il burn-in
    if (i > num_adaptive && i <= num_burn_in) {
      t <- i - num_adaptive
      theta <- curr_state$theta
      curr_eta <- c(
        log(theta$kappa),
        log(theta$sigma_beta),
        theta$v1,
        theta$v2
      )

      if (t == 1L) {
        eta_mean <- curr_eta
      } else {
        diff_eta <- curr_eta - eta_mean
        eta_cov_mat <- ((t - 2) / (t - 1)) *
          eta_cov_mat +
          (1 / t) * tcrossprod(diff_eta)
        eta_mean <- eta_mean + diff_eta / t
      }

      if (t > 10L) {
        step_rate <- opt_rate * eta_cov_mat + 1e-6 * diag(4)
      }
    }

    # 6. Saving full samples for theta and sigma2_epsilon
    theta_samples[, i] <- c(
      curr_state$theta$kappa,
      curr_state$theta$sigma_beta,
      curr_state$theta$v1,
      curr_state$theta$v2
    )
    sigma2_samples[i] <- curr_state$sigma2_epsilon

    # Saving samples after burn-in
    if (i > num_burn_in) {
      idx <- i - num_burn_in
      beta_samples[,, idx] <- curr_state$beta_mat
      mu_samples[, idx] <- curr_state$mu
      if (isTRUE(res_mh$accepted)) {
        num_accepted <- num_accepted + 1L
      }
    }

    # Updating the progress bar after every 10 iterations
    if (verbose && (i %% 10L == 0L || i == num_iterations)) {
      setTxtProgressBar(pb, i)
    }
  }

  list(
    beta_samples = beta_samples,
    mu_samples = mu_samples,
    theta_samples = theta_samples,
    sigma2_samples = sigma2_samples,
    acceptance_rate = num_accepted / num_saved,
    final_step_rate = step_rate
  )
}
