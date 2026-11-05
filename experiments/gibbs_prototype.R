# ==============================================================================
# GIBBS SAMPLER PROTOTYPE (OPTIMIZED MULTI-YEAR: i = 1, ..., M | ISOTROPO)
# ==============================================================================
devtools::load_all()

# Setting up the observations
set.seed(42)
x <- seq(0, 5, length.out = 10)
y <- seq(0, 5, length.out = 10)
locs <- as.matrix(expand.grid(x = x, y = y))


# Triangulating the domain and computing the mesh
mesh <- fm_mesh_2d(loc = locs)

N <- mesh$n
M <- 10
cat("Nr. of vertices (N):", N, "| Years (M):", M, "\n")
plot(mesh)

# Building isotropic model for the sampler
matrices <- build_spde_matrices(mesh, manifold = "planar")
K_iso <- matrices$K
C_mat <- matrices$C
diag_C <- Matrix::diag(C_mat)


# Data Generation
sim <- generate_toy_data(
  mesh = mesh,
  num_year = M,
  a_prior_mean = 5.0,
  a_prior_sd = 2.0,
  mu_prior_mean = 15.0,
  Q_mu = Diagonal(N, 1),
  sigma_e = sqrt(0.25)
)

Y_mat <- sim$y_mat
mu_true <- sim$mu_true
a_true_vec <- sim$a_true
beta_true_mat <- sim$beta_tilde_mat


# Prior specifications for the models
prior_a_prec <- 1 / 100
prior_a_mean <- 0
prior_ig_a <- 0.01
prior_ig_b <- 0.01
prior_mu_prec <- Diagonal(N, 1 / 100)
prior_mu_mean <- rep(0, N)

# MCMC Initialization
G <- 10000
samples_a1 <- numeric(G)
samples_mu1 <- numeric(G)
samples_sigma2 <- numeric(G)
samples_tau2 <- numeric(G)
samples_kappa2 <- numeric(G)

curr_a <- rep(0, M)
curr_mu <- rep(0, N)
curr_beta <- matrix(0, nrow = N, ncol = M)
curr_sigma2 <- 1.0
curr_tau2 <- 0.5
curr_kappa2 <- 0.5

# Pre-computing the SPDE Matrix
core_mat <- curr_kappa2 * C_mat + K_iso
Q_prior_beta <- forceSymmetric(
  curr_tau2 * core_mat %*% (Diagonal(x = 1 / diag_C) %*% core_mat)
)

# 5. LOOP GIBBS OTTIMIZZATO
cat("Starting Gibbs Sampler...\n")
pb <- txtProgressBar(min = 1, max = G, style = 3)


start_time <- Sys.time()
for (g in 1:G) {
  # --- STEP 1: updating a_i ---
  res_a <- update_a(
    N,
    prior_a_mean,
    prior_a_prec,
    Y_mat,
    curr_mu,
    curr_beta,
    curr_sigma2
  )

  # --- STEP 2: updating mu---
  res_mu_sum <- rowSums(Y_mat - rep(curr_a, each = N) - curr_beta)
  Q_post_mu <- prior_mu_prec + Diagonal(N, M / curr_sigma2)

  U_mu <- chol(Q_post_mu)
  rhs_mu <- as.vector(
    prior_mu_prec %*% prior_mu_mean + (res_mu_sum / curr_sigma2)
  )
  mu_post_mean <- as.vector(solve(U_mu, solve(t(U_mu), rhs_mu)))
  curr_mu <- as.vector(mu_post_mean + solve(U_mu, rnorm(N)))

  # --- STEP 3: Updating beta_tilde_i ---
  Q_post_beta <- Q_prior_beta + Diagonal(N, 1 / curr_sigma2)
  U_beta <- chol(Q_post_beta)

  RHS_beta <- (Y_mat - rep(curr_a, each = N) - curr_mu) / curr_sigma2
  Z_mat <- matrix(rnorm(N * M), nrow = N, ncol = M)

  Mean_beta <- solve(U_beta, solve(t(U_beta), RHS_beta))
  curr_beta <- as.matrix(Mean_beta + solve(U_beta, Z_mat))

  # --- STEP 4: Updating sigma2_epsilon ---
  eps_mat <- Y_mat - rep(curr_a, each = N) - curr_mu - curr_beta
  a_post_sig <- (N * M) / 2 + prior_ig_a
  b_post_sig <- sum(eps_mat^2) / 2 + prior_ig_b
  curr_sigma2 <- 1 / rgamma(1, shape = a_post_sig, rate = b_post_sig)

  # --- STEP 5: Metropolis-Withing-Gibbs ---
  log_tau2_prop <- rnorm(1, log(curr_tau2), sd = 0.08)
  log_kappa2_prop <- rnorm(1, log(curr_kappa2), sd = 0.08)
  tau2_prop <- exp(log_tau2_prop)
  kappa2_prop <- exp(log_kappa2_prop)

  core_prop <- kappa2_prop * C_mat + K_iso
  Q_prop <- forceSymmetric(
    tau2_prop * core_prop %*% (Diagonal(x = 1 / diag_C) %*% core_prop)
  )

  chol_curr <- chol(Q_prior_beta)
  chol_prop <- chol(Q_prop)

  log_det_curr <- 2 * sum(log(diag(chol_curr)))
  log_det_prop <- 2 * sum(log(diag(chol_prop)))

  quad_curr <- sum((chol_curr %*% curr_beta)^2)
  quad_prop <- sum((chol_prop %*% curr_beta)^2)

  log_lik_curr <- 0.5 * M * log_det_curr - 0.5 * quad_curr
  log_lik_prop <- 0.5 * M * log_det_prop - 0.5 * quad_prop

  log_prior_curr <- dnorm(log(curr_tau2), 0, 5, log = TRUE) +
    dnorm(log(curr_kappa2), 0, 5, log = TRUE)
  log_prior_prop <- dnorm(log_tau2_prop, 0, 5, log = TRUE) +
    dnorm(log_kappa2_prop, 0, 5, log = TRUE)

  if (
    log(runif(1)) <
      (log_lik_prop + log_prior_prop - log_lik_curr - log_prior_curr)
  ) {
    curr_tau2 <- tau2_prop
    curr_kappa2 <- kappa2_prop
    Q_prior_beta <- Q_prop
  }

  # Saving
  samples_a1[g] <- curr_a[1]
  samples_mu1[g] <- curr_mu[1]
  samples_sigma2[g] <- curr_sigma2
  samples_tau2[g] <- curr_tau2
  samples_kappa2[g] <- curr_kappa2

  setTxtProgressBar(pb, g)
}

end_time <- Sys.time()
close(pb)

elapsed_time <- end_time - start_time
cat("\nElapsed Time:\n")
print(elapsed_time)

# TRACEPLOT
par(mfrow = c(2, 2))
plot(samples_a1, type = "l", main = "a_1 (Anno 1)")
abline(h = a_true_vec[1], col = "red", lwd = 2)
plot(samples_sigma2, type = "l", main = "sigma2 (Errore)")
abline(h = 0.25, col = "red", lwd = 2)
plot(samples_tau2, type = "l", main = "tau^2")
abline(h = 1, col = "red", lwd = 2)
plot(samples_kappa2, type = "l", main = "kappa^2")
abline(h = 1, col = "red", lwd = 2)


burnin <- 1000
idx_post <- (burnin + 1):G

media_a1 <- mean(samples_a1[idx_post])
media_sigma2 <- mean(samples_sigma2[idx_post])
media_tau2 <- mean(samples_tau2[idx_post])
media_kappa2 <- mean(samples_kappa2[idx_post])

cat("\n=== CONFRONTO STIME POSTERIORI VS VERITÀ ===\n")
cat(sprintf("a_1     -> Stimato: %.4f | Vero: %.4f\n", media_a1, a_true_vec[1]))
cat(sprintf(
  "sigma^2 -> Stimato: %.4f | Vero: %.4f\n",
  media_sigma2,
  0.25
))
cat(sprintf("tau^2   -> Stimato: %.4f | Vero: %.4f\n", media_tau2, 1.0))
cat(sprintf(
  "kappa^2 -> Stimato: %.4f | Vero: %.4f\n",
  media_kappa2,
  1.0
))

cat("\nIntervalli di credibilità al 95%:\n")
cat("tau^2:", quantile(samples_tau2[idx_post], c(0.025, 0.975)), "\n")
cat("kappa^2:", quantile(samples_kappa2[idx_post], c(0.025, 0.975)), "\n")


update_a <- function(
  N,
  mean_a_old,
  tau_a_old,
  Y_mat,
  curr_mu,
  curr_beta,
  curr_sigma2
) {
  new_a <- numeric(N)
  tau_a_data <- N / curr_sigma2
  tau_a_new <- tau_a_old + tau_a_data
  for (i in 1:M) {
    mean_a_data_i <- mean(Y_mat[, i] - curr_mu - curr_beta[, i])
    mean_a_i_new <- (tau_a_old * mean_a_old[i] + tau_a_data * mean_a_data_i) /
      tau_a_new
    new_a[i] <- rnorm(1, mean = mean_a_i_new, sd = sqrt(1 / tau_a_new))
  }
  return(list(
    new_a = new_a,
    tau_a_new = tau_a_new
  ))
}
