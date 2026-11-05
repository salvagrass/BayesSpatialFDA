update_a <- function(
  N,
  M,
  mu_a_old,
  tau_a_old,
  Y_mat,
  curr_mu,
  curr_beta,
  curr_sigma2
) {
  new_samples <- numeric(M)
  tau_data <- N / curr_sigma2
  tau_a_new <- tau_a_old + tau_data

  residual_mat <- Y_mat - curr_mu - curr_beta
  mean_data <- colMeans(residual_mat)
  mean_a_new <- (tau_a_old * mu_a_old + mean_data * tau_data) / tau_a_new
  new_a <- rnorm(M, mean_a_new, sqrt(1 / tau_a_new))

  return(new_a)
}

update_mu <- function(
  N,
  M,
  mean_mu_old,
  mean_q_old,
  Y_mat,
  curr_a,
  curr_beta,
  curr_sigma2
) {
  q_data <- Diagonal(N, M / curr_sigma2)
  mean_data <- 1 / M * rowSums(Y_mat - rep(curr_a, each = N) - curr_beta)

  new_q <- q_data + prior_mu_prec
  rhs <- as.vector(q_data %*% mean_data + mean_q_old %*% mean_mu_old)

  # Speedy way to compute the inverse
  chol_factor <- chol(new_q)
  new_mean <- solve(
    chol_factor,
    solve(t(chol_factor), rhs)
  )
  new_mean <- as.vector(new_mean)
  new_mu <- sample_normal_vector(
    new_mean,
    new_q,
    is_precision = TRUE
  )
  return(new_mu)
}


update_beta <- function(
  N,
  M,
  beta_prior_prec,
  Y_mat,
  curr_a,
  curr_mu,
  curr_sigma2
) {
  new_q <- beta_prior_prec + Diagonal(N, 1 / curr_sigma2)
  rhs <- (Y_mat - rep(curr_a, each = N) - curr_mu) / curr_sigma2

  chol_factor <- chol(new_q)
  new_mean <- solve(chol_factor, solve(t(chol_factor), rhs))

  # Computing the sampling matrix, with just one Cholesky factorization
  Z_mat <- matrix(rnorm(N * M), nrow = N, ncol = M)
  curr_beta <- as.matrix(new_mean + solve(chol_factor, Z_mat))

  return(curr_beta)
}


#' Metropolis-Hastings step for SPDE hyperparameters (tau^2, kappa^2)
metropolis_hastings <- function(
  curr_tau2,
  curr_kappa2,
  curr_beta,
  Q_prior_beta,
  C_mat,
  K_iso,
  diag_C,
  prop_sd = 0.08
) {
  M <- ncol(curr_beta)

  log_tau2_prop <- rnorm(1, log(curr_tau2), sd = prop_sd)
  log_kappa2_prop <- rnorm(1, log(curr_kappa2), sd = prop_sd)
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

  # 8. Accept/Reject
  accepted <- FALSE
  if (
    log(runif(1)) <
      (log_lik_prop + log_prior_prop - log_lik_curr - log_prior_curr)
  ) {
    curr_tau2 <- tau2_prop
    curr_kappa2 <- kappa2_prop
    Q_prior_beta <- Q_prop
    accepted <- TRUE
  }

  return(list(
    tau2 = curr_tau2,
    kappa2 = curr_kappa2,
    Q_beta = Q_prior_beta,
    accepted = accepted
  ))
}


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
  a_prior_mean = 5,
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
  curr_a <- update_a(
    N,
    M,
    prior_a_mean,
    prior_a_prec,
    Y_mat,
    curr_mu,
    curr_beta,
    curr_sigma2
  )

  # centering curr_a
  a_bar <- mean(curr_a)
  curr_a <- curr_a - a_bar
  curr_mu <- curr_mu + a_bar

  # --- STEP 2: updating mu---
  curr_mu <- update_mu(
    N,
    M,
    prior_mu_mean,
    prior_mu_prec,
    Y_mat,
    curr_a,
    curr_beta,
    curr_sigma2
  )

  # --- STEP 3: Updating beta_tilde_i ---
  curr_beta <- update_beta(
    N,
    M,
    Q_prior_beta,
    Y_mat,
    curr_a,
    curr_mu,
    curr_sigma2
  )

  # Centering beta
  beta_means <- colMeans(curr_beta)
  curr_beta <- curr_beta - rep(beta_means, each = N)
  curr_a <- curr_a + beta_means

  # --- STEP 4: Updating sigma2_epsilon ---
  eps_mat <- Y_mat - rep(curr_a, each = N) - curr_mu - curr_beta
  a_post_sig <- (N * M) / 2 + prior_ig_a
  b_post_sig <- sum(eps_mat^2) / 2 + prior_ig_b
  curr_sigma2 <- 1 / rgamma(1, shape = a_post_sig, rate = b_post_sig)

  # --- STEP 5: Metropolis-Withing-Gibbs ---
  mh_step <- metropolis_hastings(
    curr_tau2 = curr_tau2,
    curr_kappa2 = curr_kappa2,
    curr_beta = curr_beta,
    Q_prior_beta = Q_prior_beta,
    C_mat = C_mat,
    K_iso = K_iso,
    diag_C = diag_C,
    prop_sd = 0.08
  )

  curr_tau2 <- mh_step$tau2
  curr_kappa2 <- mh_step$kappa2
  Q_prior_beta <- mh_step$Q_beta

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


a_true_centered <- a_true_vec - mean(a_true_vec)

cat("\n=== CONFRONTO STIME POSTERIORI VS VERITÀ ===\n")
cat(sprintf(
  "a_1     -> Stimato: %.4f | Vero: %.4f\n",
  media_a1,
  a_true_centered[1]
))
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
