update_sigma2_epsilon <- function(
  N, # Number of observations
  V, # Number of vertices
  M, # Number of years
  Phi, # Projector Matrix
  Y_mat, # N x M matrix of observations
  beta_mat, # V x M matrix of coefficients (beta in this case)
  priors # list of priors parameters
) {
  # Compute the residuals matrix and its norm
  res <- Y_mat - Phi %*% beta_mat
  norm_res <- sum(res^2)

  # Compute the new parameters for the posterior
  a_post <- priors$sigma2_epsilon$a + M * N * 0.5
  b_post <- priors$sigma2_epsilon$b + 0.5 * norm_res

  1 / rgamma(1, a_post, rate = b_post)
}


update_mu_centered <- function(
  V, # Number of vertices
  M, # Number of years
  beta_mat,
  spde_theta, # theta parameters for the precision matrix,
  spde_matrices,
  priors
) {
  # Computing the prior matrix
  Q_prior <- (1 / priors$mu$sigma2_mu) * Diagonal(V)
  Q_spde <- build_precision_matrix(
    spde_matrices,
    spde_theta$kappa,
    spde_theta$sigma_beta,
    spde_theta$v1,
    spde_theta$v2
  )
  # Computing the post parameters
  Q_post <- Q_prior + M * Q_spde
  rhs <- (1 / priors$mu$sigma2_mu) *
    priors$mu$mean +
    Q_spde %*% rowSums(beta_mat)
  mean_post <- as.vector(solve(Q_post, rhs))

  sample_normal_vector(
    mean = mean_post,
    mat = Q_post,
    N = V,
    is_precision = TRUE
  )
}

asis_step <- function(
  V,
  M,
  Phi,
  Y_mat,
  beta_mat,
  mu_old,
  sigma2_epsilon,
  priors
) {
  tilde_beta_mat <- beta_mat - mu_old

  Q_post <- (1 / priors$mu$sigma2_mu) *
    Diagonal(V) +
    M / sigma2_epsilon * t(Phi) %*% Phi

  # Computing the residuals
  res_sum <- rowSums(Y_mat - Phi %*% tilde_beta_mat)
  rhs <- (1 / priors$mu$sigma2_mu) *
    priors$mu$mean +
    (1 / sigma2_epsilon) * t(Phi) %*% res_sum

  mean_post <- as.vector(solve(Q_post, rhs))
  asis_mu <- sample_normal_vector(mean_post, Q_post, V, is_precision = TRUE)
  list(
    beta_mat = tilde_beta_mat + asis_mu,
    mu = asis_mu
  )
}


metropolis_step <- function(
  N,
  V,
  M,
  Phi,
  step_rate,
  spde_theta,
  mu,
  Y_mat,
  sigma2_epsilon,
  priors,
  spde_matrices
) {
  curr_state <- c(
    log(spde_theta$kappa),
    log(spde_theta$sigma_beta),
    spde_theta$v1,
    spde_theta$v2
  )

  new_state <- as.vector(sample_normal_vector(
    curr_state,
    step_rate
  ))

  pi_old <- .compute_log_likelihood(
    M,
    Phi,
    curr_state,
    mu,
    Y_mat,
    sigma2_epsilon,
    priors,
    spde_matrices
  )
  pi_new <- .compute_log_likelihood(
    M,
    Phi,
    new_state,
    mu,
    Y_mat,
    sigma2_epsilon,
    priors,
    spde_matrices
  )

  if (log(runif(1)) <= pi_new$value - pi_old$value) {
    new_theta <- list(
      kappa = exp(new_state[1]),
      sigma_beta = exp(new_state[2]),
      v1 = new_state[3],
      v2 = new_state[4]
    )
    is_accepted <- TRUE
    L_post <- pi_new$L_beta_post
    w_mat <- pi_new$w_mat
  } else {
    new_theta <- spde_theta
    is_accepted <- FALSE
    L_post <- pi_old$L_beta_post
    w_mat <- pi_old$w_mat
  }
  new_beta <- .update_beta(
    V,
    M,
    L_post,
    w_mat
  )
  list(
    spde_theta = new_theta,
    beta_mat = new_beta,
    accepted = is_accepted
  )
}


.compute_log_likelihood <- function(
  M,
  Phi,
  state,
  mu,
  Y_mat,
  sigma2_epsilon,
  priors,
  spde_matrices
) {
  kappa <- exp(state[1])
  sigma_beta <- exp(state[2])

  Q <- build_precision_matrix(
    spde_matrices,
    kappa,
    sigma_beta,
    state[3],
    state[4]
  )

  Q_beta_post <- Q +
    1 /
      (sigma2_epsilon) *
      t(Phi) %*% Phi

  # Compute the log determinants
  L <- chol(Q)
  L_post <- chol(Q_beta_post)
  log_det_Q <- 2 * sum(log(diag(L)))
  log_det_Q_post <- 2 * sum(log(diag(L_post)))

  # Compute the contributions of the quadratic forms
  quad_form_Q <- sum((L %*% mu)^2)

  linear_mat <- as.vector(Q %*% mu) + (1 / sigma2_epsilon) * t(Phi) %*% Y_mat
  w_mat <- solve(t(L_post), linear_mat)
  quad_form_post <- sum((w_mat)^2)

  # Compute the contributions from the priors
  prior_fact <- -priors$lambda_kappa *
    kappa -
    priors$lambda_sigma_beta * sigma_beta -
    0.5 * (state[3]^2 + state[4]^2) / priors$sigma_v^2

  # Compute the logarithm of the jacobian determinant
  scale_fact <- state[1] + state[2]

  # Compute the total likelihood
  pi <- prior_fact +
    scale_fact +
    M * 0.5 * log_det_Q -
    M * 0.5 * log_det_Q_post -
    M * 0.5 * quad_form_Q +
    0.5 * quad_form_post

  list(
    value = pi,
    Q_beta_post = Q_beta_post,
    L_beta_post = L_post,
    w_mat = w_mat
  )
}

.update_beta <- function(
  V,
  M,
  L_post,
  w_mat
) {
  # w_mat = L_post^(-T)*linear_mat, mean_mat = L_post^(-1)*L_post^(-T)*linear_mat
  # so mean_mat = L_post^(-1)*w_mat
  Z <- matrix(rnorm(V * M), nrow = V, ncol = M)
  as.matrix(solve(L_post, w_mat + Z))
}
