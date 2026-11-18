# experiments/gibbs_sampler/mcmc_updates.R
# Core MCMC update functions for the Marginalised Anisotropic SPDE Gibbs Sampler

#' Build the Anisotropic SPDE Prior Precision Matrix Q_pr
#'
#' @param tau2 Precision scaling parameter
#' @param kappa2 Spatial range scaling parameter (kappa^2)
#' @param gamma Anisotropy ratio (> 0)
#' @param theta Anisotropy rotation angle in [0, pi)
#' @param C_mat Mass matrix C
#' @param K11,K22,K12 Stiffness matrices
#' @param inv_diag_C Vector of 1 / diag(C_mat) for mass lumping
#' @return Symmetric sparse precision matrix Q_pr
build_spde_precision <- function(
  tau2,
  kappa2,
  gamma,
  theta,
  C_mat,
  K11,
  K22,
  K12,
  inv_diag_C
) {
  cos_t <- cos(theta)
  sin_t <- sin(theta)
  inv_g <- 1 / gamma

  H11 <- gamma * (cos_t^2) + inv_g * (sin_t^2)
  H22 <- gamma * (sin_t^2) + inv_g * (cos_t^2)
  H12 <- (gamma - inv_g) * sin_t * cos_t

  G_H <- H11 * K11 + H12 * K12 + H22 * K22
  core <- kappa2 * C_mat + G_H

  Q_pr <- Matrix::forceSymmetric(
    tau2 * core %*% (Matrix::Diagonal(x = inv_diag_C) %*% core)
  )
  return(Q_pr)
}

#' Gibbs Update for the Long-Term Spatial Mean Field (mu)
#'
#' @param N Number of mesh vertices
#' @param M Number of years
#' @param prior_mean Prior mean vector for mu
#' @param prior_prec Prior precision matrix for mu
#' @param Y_mat Observation matrix (N x M)
#' @param curr_beta Current latent anomaly matrix (N x M)
#' @param curr_sigma2 Current observation noise variance
#' @return Updated mu vector of length N
update_mu <- function(
  N,
  M,
  prior_mean,
  prior_prec,
  Y_mat,
  curr_beta,
  curr_sigma2
) {
  q_data <- Matrix::Diagonal(N, M / curr_sigma2)
  mean_data <- (1 / M) * rowSums(Y_mat - curr_beta)

  new_q <- q_data + prior_prec
  rhs <- as.vector(q_data %*% mean_data + prior_prec %*% prior_mean)

  chol_factor <- Matrix::chol(new_q)
  new_mean <- as.vector(Matrix::solve(
    chol_factor,
    Matrix::solve(t(chol_factor), rhs)
  ))

  new_mu <- sample_normal_vector(new_mean, new_q, is_precision = TRUE)
  return(new_mu)
}

#' Marginalised Random-Walk Metropolis-Hastings for SPDE Hyperparameters
#'
#' Integrates out beta analytically to evaluate p(Y - mu | tau2, kappa2, gamma, theta, sigma2).
#'
#' @return List with updated hyperparameters, cached matrices/Cholesky factors, and acceptance flag
marginal_metropolis_hastings <- function(
  curr_tau2,
  curr_kappa2,
  curr_gamma,
  curr_theta,
  B_mat,
  curr_sigma2,
  C_mat,
  K11,
  K22,
  K12,
  inv_diag_C,
  prop_cov,
  Q_pr_curr,
  chol_pr_curr,
  chol_post_curr,
  M_post_curr,
  log_lik_curr,
  prior_sd = 5
) {
  M_years <- ncol(B_mat)
  N_nodes <- nrow(B_mat)

  # 1. Propose in unconstrained space
  curr_vec <- c(log(curr_tau2), log(curr_kappa2), log(curr_gamma), curr_theta)
  prop_vec <- sample_normal_vector(
    mean = curr_vec,
    mat = prop_cov,
    N = 4,
    is_precision = FALSE
  )

  tau2_prop <- exp(prop_vec[1])
  kappa2_prop <- exp(prop_vec[2])
  gamma_prop <- exp(prop_vec[3])
  theta_prop <- prop_vec[4] %% pi

  # 2. Build proposed SPDE precision and posterior precision
  Q_pr_prop <- build_spde_precision(
    tau2 = tau2_prop,
    kappa2 = kappa2_prop,
    gamma = gamma_prop,
    theta = theta_prop,
    C_mat = C_mat,
    K11 = K11,
    K22 = K22,
    K12 = K12,
    inv_diag_C = inv_diag_C
  )
  Q_post_prop <- Q_pr_prop + Matrix::Diagonal(N_nodes, 1 / curr_sigma2)

  # 3. Cholesky factorizations and posterior conditional mean M_post
  chol_pr_prop <- Matrix::chol(Q_pr_prop)
  chol_post_prop <- Matrix::chol(Q_post_prop)
  M_post_prop <- Matrix::solve(
    chol_post_prop,
    Matrix::solve(t(chol_post_prop), B_mat)
  )

  # 4. Marginal log-likelihood (with beta integrated out)
  log_det_pr_prop <- 2 * sum(log(Matrix::diag(chol_pr_prop)))
  log_det_post_prop <- 2 * sum(log(Matrix::diag(chol_post_prop)))
  quad_prop <- sum(B_mat * M_post_prop)

  log_lik_prop <- 0.5 *
    M_years *
    log_det_pr_prop -
    0.5 * M_years * log_det_post_prop +
    0.5 * quad_prop

  # 5. Log-priors on log-parameters (flat on theta in [0, pi))
  log_prior_curr <- dnorm(curr_vec[1], 0, prior_sd, log = TRUE) +
    dnorm(curr_vec[2], 0, prior_sd, log = TRUE) +
    dnorm(curr_vec[3], 0, prior_sd, log = TRUE)

  log_prior_prop <- dnorm(prop_vec[1], 0, prior_sd, log = TRUE) +
    dnorm(prop_vec[2], 0, prior_sd, log = TRUE) +
    dnorm(prop_vec[3], 0, prior_sd, log = TRUE)

  # 6. Accept / Reject
  accepted <- FALSE
  log_alpha <- (log_lik_prop + log_prior_prop) - (log_lik_curr + log_prior_curr)

  if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
    curr_tau2 <- tau2_prop
    curr_kappa2 <- kappa2_prop
    curr_gamma <- gamma_prop
    curr_theta <- theta_prop

    Q_pr_curr <- Q_pr_prop
    chol_pr_curr <- chol_pr_prop
    chol_post_curr <- chol_post_prop
    M_post_curr <- M_post_prop
    log_lik_curr <- log_lik_prop
    accepted <- TRUE
  }

  return(list(
    tau2 = curr_tau2,
    kappa2 = curr_kappa2,
    gamma = curr_gamma,
    theta = curr_theta,
    Q_pr = Q_pr_curr,
    chol_pr = chol_pr_curr,
    chol_post = chol_post_curr,
    M_post = M_post_curr,
    log_lik = log_lik_curr,
    accepted = accepted
  ))
}
