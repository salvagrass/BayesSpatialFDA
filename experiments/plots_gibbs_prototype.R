# ==============================================================================
# GIBBS SAMPLER PROTOTYPE (MARGINALIZED BLOCK-UPDATE | ADAPTIVE | ANISOTROPIC)
# Model: Y_t = mu + beta_t + epsilon_t
# ==============================================================================
devtools::load_all()

# ==============================================================================
# 1. GIBBS AND MARGINALIZED METROPOLIS UPDATES
# ==============================================================================

update_mu <- function(
  N,
  M,
  mean_mu_old,
  mean_q_old,
  Y_mat,
  curr_beta,
  curr_sigma2
) {
  q_data <- Diagonal(N, M / curr_sigma2)
  # Data contribution to mu is the residual after removing beta
  mean_data <- 1 / M * rowSums(Y_mat - curr_beta)

  new_q <- q_data + mean_q_old
  rhs <- as.vector(q_data %*% mean_data + mean_q_old %*% mean_mu_old)

  chol_factor <- chol(new_q)
  new_mean <- as.vector(solve(chol_factor, solve(t(chol_factor), rhs)))

  new_mu <- sample_normal_vector(new_mean, new_q, is_precision = TRUE)
  return(new_mu)
}

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
  diag_C,
  prop_cov,
  Q_pr_curr,
  chol_pr_curr,
  chol_post_curr,
  M_post_curr,
  log_lik_curr
) {
  M_years <- ncol(B_mat)
  N_nodes <- nrow(B_mat)

  # 1. Proposal in unconstrained space
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

  # 2. Build Proposed SPDE Matrices
  H11 <- gamma_prop *
    (cos(theta_prop)^2) +
    (1 / gamma_prop) * (sin(theta_prop)^2)
  H22 <- gamma_prop *
    (sin(theta_prop)^2) +
    (1 / gamma_prop) * (cos(theta_prop)^2)
  H12 <- (gamma_prop - 1 / gamma_prop) * sin(theta_prop) * cos(theta_prop)

  G_H <- H11 * K11 + H12 * K12 + H22 * K22
  core_prop <- kappa2_prop * C_mat + G_H
  Q_pr_prop <- forceSymmetric(
    tau2_prop * core_prop %*% (Diagonal(x = 1 / diag_C) %*% core_prop)
  )
  Q_post_prop <- Q_pr_prop + Diagonal(N_nodes, 1 / curr_sigma2)

  # 3. Factorizations
  chol_pr_prop <- chol(Q_pr_prop)
  chol_post_prop <- chol(Q_post_prop)
  M_post_prop <- solve(chol_post_prop, solve(t(chol_post_prop), B_mat))

  # 4. MARGINAL Log-Likelihood
  log_det_pr_prop <- 2 * sum(log(diag(chol_pr_prop)))
  log_det_post_prop <- 2 * sum(log(diag(chol_post_prop)))
  quad_prop <- sum(B_mat * M_post_prop)

  log_lik_prop <- 0.5 *
    M_years *
    log_det_pr_prop -
    0.5 * M_years * log_det_post_prop +
    0.5 * quad_prop

  # 5. Log-Priors
  log_prior_curr <- dnorm(log(curr_tau2), 0, 5, log = TRUE) +
    dnorm(log(curr_kappa2), 0, 5, log = TRUE) +
    dnorm(log(curr_gamma), 0, 5, log = TRUE)
  log_prior_prop <- dnorm(prop_vec[1], 0, 5, log = TRUE) +
    dnorm(prop_vec[2], 0, 5, log = TRUE) +
    dnorm(prop_vec[3], 0, 5, log = TRUE)

  # 6. Accept / Reject
  accepted <- FALSE
  if (
    log(runif(1)) <
      (log_lik_prop + log_prior_prop - log_lik_curr - log_prior_curr)
  ) {
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

# ==============================================================================
# 2. SETUP AND DATA GENERATION
# ==============================================================================

set.seed(42)
x <- seq(0, 5, length.out = 20)
y <- seq(0, 5, length.out = 20)
locs <- as.matrix(expand.grid(x = x, y = y))

mesh <- fm_mesh_2d(loc = locs)
N <- mesh$n
M <- 10
cat("Nr. of vertices (N):", N, "| Years (M):", M, "\n")

matrices <- build_spde_matrices(mesh, manifold = "planar")
K11 <- matrices$K11
K22 <- matrices$K22
K12 <- matrices$K12
C_mat <- matrices$C
diag_C <- Matrix::diag(C_mat)

true_tau <- 1.0
true_kappa <- 1.0
true_gamma <- 2.0
true_theta <- pi / 4
true_sigma_e <- sqrt(0.25)

# Generatore dati (assume Y_mat = beta + epsilon)
sim <- generate_toy_data(
  mesh = mesh,
  num_year = M,
  mu_prior_mean = 0.5,
  Q_mu = Diagonal(N, 1),
  tau = true_tau,
  k = true_kappa,
  gamma = true_gamma,
  theta = true_theta,
  sigma_e = true_sigma_e
)

Y_mat <- sim$y_mat
mu_true <- sim$mu_true
# L'anomalia vera spaziale (a media zero)
true_anomaly_mat <- sim$beta_mat - sim$mu_true

# ==============================================================================
# 3. MCMC INITIALIZATION
# ==============================================================================

G <- 5000
burnin <- 1000
idx_post <- (burnin + 1):G
samples_sigma2 <- numeric(G)
samples_tau2 <- numeric(G)
samples_kappa2 <- numeric(G)
samples_gamma <- numeric(G)
samples_theta <- numeric(G)
samples_beta1 <- numeric(G)

curr_mu <- rep(0, N)
curr_beta <- matrix(0, nrow = N, ncol = M) # Zero-mean starting field
curr_sigma2 <- 1.0
curr_tau2 <- 0.5
curr_kappa2 <- 0.5
curr_gamma <- 1.0
curr_theta <- 0.0

prior_ig_a <- 0.01
prior_ig_b <- 0.01
prior_mu_prec <- Diagonal(N, 1 / 100)
prior_mu_mean <- rep(0, N)

H11_init <- curr_gamma *
  (cos(curr_theta)^2) +
  (1 / curr_gamma) * (sin(curr_theta)^2)
H22_init <- curr_gamma *
  (sin(curr_theta)^2) +
  (1 / curr_gamma) * (cos(curr_theta)^2)
H12_init <- (curr_gamma - 1 / curr_gamma) * sin(curr_theta) * cos(curr_theta)
G_H_init <- H11_init * K11 + H12_init * K12 + H22_init * K22

core_mat <- curr_kappa2 * C_mat + G_H_init
Q_pr_curr <- forceSymmetric(
  curr_tau2 * core_mat %*% (Diagonal(x = 1 / diag_C) %*% core_mat)
)
chol_pr_curr <- chol(Q_pr_curr)

mcmc_history <- matrix(0, nrow = G, ncol = 4)
colnames(mcmc_history) <- c("log_tau2", "log_kappa2", "log_gamma", "theta")
prop_cov <- diag(x = 0.05^2, nrow = 4)
scale_factor <- (2.38^2) / 4 * 0.7
acc_count <- 0


beta_sum_post <- matrix(0, nrow = N, ncol = M)

# ==============================================================================
# 4. MAIN GIBBS LOOP
# ==============================================================================

cat("Starting Marginalized Adaptive Gibbs Sampler (mu + beta)...\n")
pb <- txtProgressBar(min = 1, max = G, style = 3)

start_time <- Sys.time()
for (g in 1:G) {
  # Step 1: Update long-term mean mu
  curr_mu <- update_mu(
    N,
    M,
    prior_mu_mean,
    prior_mu_prec,
    Y_mat,
    curr_beta,
    curr_sigma2
  )
  # Step 2: Update sigma2
  eps_mat <- Y_mat - curr_mu - curr_beta
  a_post_sig <- (N * M) / 2 + prior_ig_a
  b_post_sig <- sum(eps_mat^2) / 2 + prior_ig_b
  curr_sigma2 <- 1 / rgamma(1, shape = a_post_sig, rate = b_post_sig)

  # --- BLOCK UPDATE PREPARATION ---
  Y_star <- Y_mat - curr_mu # SPDE explains the residuals
  B_mat <- Y_star / curr_sigma2
  Q_post_curr <- Q_pr_curr + Diagonal(N, 1 / curr_sigma2)
  chol_post_curr <- chol(Q_post_curr)
  M_post_curr <- solve(chol_post_curr, solve(t(chol_post_curr), B_mat))

  log_det_pr_curr <- 2 * sum(log(diag(chol_pr_curr)))
  log_det_post_curr <- 2 * sum(log(diag(chol_post_curr)))
  quad_curr <- sum(B_mat * M_post_curr)
  log_lik_curr <- 0.5 *
    M *
    log_det_pr_curr -
    0.5 * M * log_det_post_curr +
    0.5 * quad_curr

  # Adaptive Logic
  if (g > 200 && (g %% 10 == 0)) {
    emp_cov <- cov(mcmc_history[1:(g - 1), ])
    prop_cov <- scale_factor * emp_cov + diag(1e-4, nrow = 4)
  }

  # Step 3: Marginal M-H
  mh_step <- marginal_metropolis_hastings(
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
    diag_C,
    prop_cov,
    Q_pr_curr,
    chol_pr_curr,
    chol_post_curr,
    M_post_curr,
    log_lik_curr
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
    acc_count <- acc_count + 1
  }

  # Step 4: Joint Beta update (Anomaly)
  Z_mat <- matrix(rnorm(N * M), nrow = N, ncol = M)
  curr_beta <- as.matrix(M_post_curr + solve(chol_post_curr, Z_mat))

  # Save parameters
  samples_sigma2[g] <- curr_sigma2
  samples_tau2[g] <- curr_tau2
  samples_kappa2[g] <- curr_kappa2
  samples_gamma[g] <- curr_gamma
  samples_theta[g] <- curr_theta
  samples_beta1[g] <- curr_beta[1, 1]
  if (g > burnin) {
    beta_sum_post <- curr_beta + beta_sum_post
  }

  mcmc_history[g, ] <- c(
    log(curr_tau2),
    log(curr_kappa2),
    log(curr_gamma),
    curr_theta
  )
  setTxtProgressBar(pb, g)
}
end_time <- Sys.time()
close(pb)

cat(sprintf(
  "\nElapsed Time: %.2f seconds\n",
  as.numeric(difftime(end_time, start_time, units = "secs"))
))
cat(sprintf("Acceptance Rate M-H: %.2f%%\n", (acc_count / G) * 100))

# ==============================================================================
# 5. DIAGNOSTICS & PLOTS
# ==============================================================================
par(mfrow = c(2, 3))
plot(samples_sigma2, type = "l", main = "sigma^2 (Noise)", xlab = "Iteration")
abline(h = true_sigma_e^2, col = "red", lwd = 2)

plot(samples_tau2, type = "l", main = "tau^2", xlab = "Iteration")
abline(h = true_tau^2, col = "red", lwd = 2)

plot(samples_kappa2, type = "l", main = "kappa^2", xlab = "Iteration")
abline(h = true_kappa^2, col = "red", lwd = 2)

plot(
  samples_gamma,
  type = "l",
  main = "gamma (Anisotropy Ratio)",
  xlab = "Iteration"
)
abline(h = true_gamma, col = "red", lwd = 2)

plot(samples_theta, type = "l", main = "theta (Angle)", xlab = "Iteration")
abline(h = true_theta, col = "red", lwd = 2)

plot(
  samples_beta1,
  type = "l",
  main = expression(beta["1,1"] ~ "Anomaly"),
  xlab = "Iteration"
)
abline(h = true_anomaly_mat[1, 1], col = "red", lwd = 2)

# ==============================================================================
# Autocorrelation plots
# ==============================================================================

post_draws <- cbind(
  "log sigma^2" = log(samples_sigma2[idx_post]),
  "log tau^2"   = log(samples_tau2[idx_post]),
  "log kappa^2" = log(samples_kappa2[idx_post]),
  "log gamma"   = log(samples_gamma[idx_post]),
  "theta"       = samples_theta[idx_post],
  "beta[1,1]"   = samples_beta1[idx_post]
)

par(mfrow = c(2, 3))
for (j in seq_len(ncol(post_draws))) {
  acf(
    post_draws[, j],
    lag.max = 60,
    main = paste("ACF:", colnames(post_draws)[j])
  )
}
par(mfrow = c(1, 1))


# ==============================================================================
# POSTERIOR ESTIMATES VS GROUND TRUTH
# ==============================================================================

cat("\n=== POSTERIOR ESTIMATES VS GROUND TRUTH ===\n")
cat(sprintf(
  "sigma^2 -> Estimated: %.4f | True: %.4f\n",
  mean(samples_sigma2[idx_post]),
  true_sigma_e^2
))
cat(sprintf(
  "tau^2   -> Estimated: %.4f | True: %.4f\n",
  mean(samples_tau2[idx_post]),
  true_tau^2
))
cat(sprintf(
  "kappa^2 -> Estimated: %.4f | True: %.4f\n",
  mean(samples_kappa2[idx_post]),
  true_kappa^2
))
cat(sprintf(
  "gamma   -> Estimated: %.4f | True: %.4f\n",
  mean(samples_gamma[idx_post]),
  true_gamma
))
cat(sprintf(
  "theta   -> Estimated: %.4f | True: %.4f\n",
  mean(samples_theta[idx_post]),
  true_theta
))
cat(sprintf(
  "beta_11 -> Estimated: %.4f | True: %.4f\n",
  mean(samples_beta1[idx_post]),
  true_anomaly_mat[1, 1]
))


beta_est_mean <- beta_sum_post / (G - burnin)

# ==============================================================================
# PLOT CONTINUO CON LA MESH E BASI FEM (Grado 1)
# ==============================================================================
par(mfrow = c(1, 2))
plot_spde_field(
  mesh,
  true_anomaly_mat[, 1],
  main_title = "True Anomaly (Year 1)"
)
plot_spde_field(
  mesh,
  beta_est_mean[, 1],
  main_title = "Estimated Anomaly (Year 1)"
)
par(mfrow = c(1, 1))
