# ==============================================================================
# GIBBS + MH SAMPLER ON REAL PLANAR CLIMATE DATA (ISOTROPIC SPDE + PC PRIORS)
# Model: Y_i = a_i * 1_N + Phi %*% (mu + beta_tilde_i) + eps_i,  i = 1, ..., M
# ==============================================================================
devtools::load_all()
library(Matrix)

# ------------------------------------------------------------------------------
# 1. HELPER FUNCTIONS
# ------------------------------------------------------------------------------

#' Build isotropic SPDE precision matrix Q(rho, sigma_beta)
build_Q_spde <- function(rho, sigma_beta, C_mat, K_iso, inv_diag_C) {
  kappa <- sqrt(8) / rho
  kappa2 <- kappa^2
  tau2 <- 1 / (4 * pi * kappa2 * sigma_beta^2)
  core <- kappa2 * C_mat + K_iso
  forceSymmetric(tau2 * (core %*% (inv_diag_C %*% core)))
}

#' Sample x ~ N(Q^{-1} rhs, Q^{-1}) using ONE AMD-permuted sparse Cholesky
#' Handles both vector rhs (for mu) and V x M matrix rhs (for beta)
sample_gmrf_ch <- function(ch, rhs) {
  mean_part <- solve(ch, rhs, system = "A")
  if (is.matrix(rhs) || inherits(rhs, "Matrix")) {
    Z <- matrix(rnorm(nrow(rhs) * ncol(rhs)), nrow = nrow(rhs), ncol = ncol(rhs))
    y <- solve(ch, Z, system = "Lt")
    noise <- solve(ch, y, system = "Pt")
    return(as.matrix(mean_part + noise))
  } else {
    Z <- rnorm(length(rhs))
    y <- solve(ch, Z, system = "Lt")
    noise <- solve(ch, y, system = "Pt")
    return(as.vector(mean_part + noise))
  }
}

#' Log-determinant log|Q| from a CHMfactor (P Q P^T = L L^T)
log_det_ch <- function(ch) {
  2 * as.numeric(determinant(ch, logarithm = TRUE, sqrt = TRUE)$modulus)
}

#' Step 1: Full conditional for global yearly effects a_1, ..., a_M
update_a <- function(N, M, mu_a_prior, tau_a_prior, Y_mat, Phi, curr_mu, curr_beta, curr_sigma2) {
  Phi_field <- as.matrix(Phi %*% (curr_beta + curr_mu))
  residual_mat <- Y_mat - Phi_field
  mean_data <- colMeans(residual_mat)

  tau_data <- N / curr_sigma2
  tau_a_post <- tau_a_prior + tau_data
  mean_a_post <- (tau_a_prior * mu_a_prior + tau_data * mean_data) / tau_a_post

  rnorm(M, mean_a_post, sqrt(1 / tau_a_post))
}

#' Step 2: Full conditional for spatial mean field mu (length V)
update_mu <- function(M, prior_mu_mean, prior_mu_prec, Y_mat, Phi, PhiT, PhiTPhi,
                      curr_a, curr_beta, curr_sigma2) {
  # Sum of residuals over all M years: sum_i (Y_i - a_i - Phi %*% beta_i)
  res_sum <- rowSums(Y_mat) - sum(curr_a) - as.vector(Phi %*% rowSums(curr_beta))

  Q_post <- forceSymmetric(prior_mu_prec + (M / curr_sigma2) * PhiTPhi)
  rhs <- as.vector(prior_mu_prec %*% prior_mu_mean + (1 / curr_sigma2) * (PhiT %*% res_sum))

  ch <- Cholesky(Q_post, LDL = FALSE, perm = TRUE)
  sample_gmrf_ch(ch, rhs)
}

#' Step 3: Full conditional for yearly spatial anomalies beta_tilde (V x M)
update_beta <- function(N, Q_prior_beta, Y_mat, Phi, PhiT, PhiTPhi,
                        curr_a, curr_mu, curr_sigma2) {
  Phi_mu <- as.vector(Phi %*% curr_mu)
  res_mat <- Y_mat - rep(curr_a, each = N) - Phi_mu

  Q_post <- forceSymmetric(Q_prior_beta + (1 / curr_sigma2) * PhiTPhi)
  rhs <- as.matrix((1 / curr_sigma2) * (PhiT %*% res_mat))

  ch <- Cholesky(Q_post, LDL = FALSE, perm = TRUE)
  sample_gmrf_ch(ch, rhs)
}

#' PC log-prior density in transformed (log_rho, log_sigma) coordinates
#' Includes Jacobian (+log_rho + log_sigma) for symmetric RW proposal on log scale
log_pc_prior_eta <- function(log_rho, log_sigma, lambda_rho, lambda_sigma) {
  rho <- exp(log_rho)
  sigma <- exp(log_sigma)
  -log_rho - lambda_rho / rho + log_sigma - lambda_sigma * sigma
}

#' Step 5: Metropolis-Hastings for (log_rho, log_sigma_beta) with PC priors
metropolis_hastings_pc <- function(
    curr_rho,
    curr_sigma_beta,
    curr_beta,
    Q_prior_beta,
    ch_prior_beta,
    C_mat,
    K_iso,
    inv_diag_C,
    lambda_rho,
    lambda_sigma,
    prop_sd = c(0.12, 0.12)
) {
  M <- ncol(curr_beta)

  log_rho_curr <- log(curr_rho)
  log_sig_curr <- log(curr_sigma_beta)

  log_rho_prop <- rnorm(1, log_rho_curr, sd = prop_sd[1])
  log_sig_prop <- rnorm(1, log_sig_curr, sd = prop_sd[2])
  rho_prop <- exp(log_rho_prop)
  sig_prop <- exp(log_sig_prop)

  # Guard against extreme numerical proposals
  if (rho_prop < 1e-4 || rho_prop > 1e4 || sig_prop < 1e-4 || sig_prop > 1e3) {
    return(list(
      rho = curr_rho, sigma_beta = curr_sigma_beta,
      Q_beta = Q_prior_beta, ch_beta = ch_prior_beta, accepted = FALSE
    ))
  }

  Q_prop <- build_Q_spde(rho_prop, sig_prop, C_mat, K_iso, inv_diag_C)
  ch_prop <- tryCatch(
    Cholesky(Q_prop, LDL = FALSE, perm = TRUE),
    error = function(e) NULL
  )
  if (is.null(ch_prop)) {
    return(list(
      rho = curr_rho, sigma_beta = curr_sigma_beta,
      Q_beta = Q_prior_beta, ch_beta = ch_prior_beta, accepted = FALSE
    ))
  }

  log_det_curr <- log_det_ch(ch_prior_beta)
  log_det_prop <- log_det_ch(ch_prop)

  quad_curr <- sum(curr_beta * as.matrix(Q_prior_beta %*% curr_beta))
  quad_prop <- sum(curr_beta * as.matrix(Q_prop %*% curr_beta))

  log_lik_curr <- 0.5 * M * log_det_curr - 0.5 * quad_curr
  log_lik_prop <- 0.5 * M * log_det_prop - 0.5 * quad_prop

  log_prior_curr <- log_pc_prior_eta(log_rho_curr, log_sig_curr, lambda_rho, lambda_sigma)
  log_prior_prop <- log_pc_prior_eta(log_rho_prop, log_sig_prop, lambda_rho, lambda_sigma)

  log_alpha <- (log_lik_prop + log_prior_prop) - (log_lik_curr + log_prior_curr)

  accepted <- FALSE
  if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
    curr_rho <- rho_prop
    curr_sigma_beta <- sig_prop
    Q_prior_beta <- Q_prop
    ch_prior_beta <- ch_prop
    accepted <- TRUE
  }

  list(
    rho = curr_rho,
    sigma_beta = curr_sigma_beta,
    Q_beta = Q_prior_beta,
    ch_beta = ch_prior_beta,
    accepted = accepted
  )
}

# ------------------------------------------------------------------------------
# 2. PREPARE REAL DATA (2000 x 2000 km SQUARE CENTERED ON ROME, UNIT-SCALED)
# ------------------------------------------------------------------------------
set.seed(42)

real_data <- prepare_real_planar_data(
  path = "data/climate_data.RData",
  center = c(41.9028, 12.4964),
  shape = "square",
  side_km = 2000,
  scale_to_unit = TRUE,
  build_mesh = TRUE,
  max_edge = c(0.08, 0.25),
  offset = c(0.02, 0.25) # Thin inner rim + coarse outer buffer
)

Y_mat <- real_data$y_mat
Phi   <- real_data$Phi
mesh  <- real_data$mesh
years <- real_data$years
scale_km <- real_data$meta$scale_factor_km

N <- nrow(Y_mat)
M <- ncol(Y_mat)
V <- mesh$n

cat(sprintf("Real data loaded -> N obs: %d | V mesh nodes: %d | M years: %d (%d-%d)\n",
            N, V, M, min(years), max(years)))

# Precompute FEM matrices and cross-products ONCE
spde_mats  <- build_spde_matrices(mesh, manifold = "planar")
C_mat      <- spde_mats$C
K_iso      <- spde_mats$K
diag_C     <- Matrix::diag(C_mat)
inv_diag_C <- Diagonal(x = 1 / diag_C)

PhiT    <- t(Phi)
PhiTPhi <- forceSymmetric(crossprod(Phi))

# ------------------------------------------------------------------------------
# 3. PRIORS SPECIFICATION
# ------------------------------------------------------------------------------
# Global yearly effects a_i ~ N(0, 10^2)
prior_a_mean <- 0.0
prior_a_prec <- 1 / 100

# Observation error variance sigma_eps^2 ~ Inv-Gamma(0.01, 0.01)
prior_ig_a <- 0.01
prior_ig_b <- 0.01

# Spatial prior on long-term mean field mu ~ N_V(mu_0, Q_mu^{-1})
# Weakly informative Matérn prior: range ~ 0.5 (1000 km), marginal SD ~ 20 degC
prior_mu_mean <- rep(mean(Y_mat), V)
prior_mu_prec <- build_Q_spde(rho = 0.5, sigma_beta = 20.0, C_mat, K_iso, inv_diag_C)

# PC Priors on SPDE anomalies (rho in unit scale [1 unit = 2000 km], sigma_beta in degC):
#   P(rho < 0.10) = 0.05  -> 5% prior prob range < 200 km
#   P(sigma_beta > 5.0) = 0.05 -> 5% prior prob yearly spatial anomaly SD > 5 degC
rho0 <- 0.10; alpha_rho <- 0.05
sig0 <- 5.0;  alpha_sig <- 0.05
lambda_rho   <- -rho0 * log(alpha_rho)
lambda_sigma <- -log(alpha_sig) / sig0

# ------------------------------------------------------------------------------
# 4. MCMC INITIALIZATION & STORAGE
# ------------------------------------------------------------------------------
G      <- 5000L
burnin <- 1500L
n_post <- G - burnin

samples_a          <- matrix(0, nrow = G, ncol = M)
samples_sigma2     <- numeric(G)
samples_rho        <- numeric(G)
samples_sigma_beta <- numeric(G)
mh_accepted        <- logical(G)

# Running posterior summaries (after burn-in) for full spatial fields
mu_sum    <- numeric(V)
mu_sq_sum <- numeric(V)
# Save reconstructed field for one target year (e.g., 2003 heatwave or 2007)
target_year_idx <- which(years == 2007L)
f_target_sum    <- numeric(N)
f_target_sq_sum <- numeric(N)

# Initial values
curr_a          <- colMeans(Y_mat) - mean(Y_mat)
curr_mu         <- rep(mean(Y_mat), V)
curr_beta       <- matrix(0, nrow = V, ncol = M)
curr_sigma2     <- 0.5
curr_rho        <- 0.5   # 1000 km
curr_sigma_beta <- 1.0   # 1 degC

Q_prior_beta  <- build_Q_spde(curr_rho, curr_sigma_beta, C_mat, K_iso, inv_diag_C)
ch_prior_beta <- Cholesky(Q_prior_beta, LDL = FALSE, perm = TRUE)

# ------------------------------------------------------------------------------
# 5. GIBBS + METROPOLIS-HASTINGS LOOP
# ------------------------------------------------------------------------------
cat("Starting Gibbs + MH Sampler on real data...\n")
pb <- txtProgressBar(min = 1, max = G, style = 3)
start_time <- Sys.time()

for (g in seq_len(G)) {
  # --- STEP 1: Update a_i ---
  curr_a <- update_a(
    N, M, prior_a_mean, prior_a_prec,
    Y_mat, Phi, curr_mu, curr_beta, curr_sigma2
  )

  # Center (a, mu): leaves a_i * 1_N + Phi %*% mu strictly unchanged
  a_bar <- mean(curr_a)
  curr_a  <- curr_a - a_bar
  curr_mu <- curr_mu + a_bar

  # --- STEP 2: Update mu ---
  curr_mu <- update_mu(
    M, prior_mu_mean, prior_mu_prec,
    Y_mat, Phi, PhiT, PhiTPhi,
    curr_a, curr_beta, curr_sigma2
  )

  # --- STEP 3: Update beta_tilde_i ---
  curr_beta <- update_beta(
    N, Q_prior_beta,
    Y_mat, Phi, PhiT, PhiTPhi,
    curr_a, curr_mu, curr_sigma2
  )

  # --- STEP 4: Update sigma2_epsilon ---
  Phi_total <- as.matrix(Phi %*% (curr_beta + curr_mu))
  eps_mat   <- Y_mat - rep(curr_a, each = N) - Phi_total
  a_post_sig <- (N * M) / 2 + prior_ig_a
  b_post_sig <- sum(eps_mat^2) / 2 + prior_ig_b
  curr_sigma2 <- 1 / rgamma(1, shape = a_post_sig, rate = b_post_sig)

  # --- STEP 5: MH step for (rho, sigma_beta) BEFORE recentering beta ---
  mh_step <- metropolis_hastings_pc(
    curr_rho        = curr_rho,
    curr_sigma_beta = curr_sigma_beta,
    curr_beta       = curr_beta,
    Q_prior_beta    = Q_prior_beta,
    ch_prior_beta   = ch_prior_beta,
    C_mat           = C_mat,
    K_iso           = K_iso,
    inv_diag_C      = inv_diag_C,
    lambda_rho      = lambda_rho,
    lambda_sigma    = lambda_sigma,
    prop_sd         = c(0.10, 0.10)
  )
  curr_rho        <- mh_step$rho
  curr_sigma_beta <- mh_step$sigma_beta
  Q_prior_beta    <- mh_step$Q_beta
  ch_prior_beta   <- mh_step$ch_beta
  mh_accepted[g]  <- mh_step$accepted

  # Recenter (beta, a) for identifiability & mixing: leaves a_i + Phi %*% beta_i unchanged
  beta_means <- colMeans(curr_beta)
  curr_beta  <- curr_beta - rep(beta_means, each = V)
  curr_a     <- curr_a + beta_means

  # Re-enforce mean(a) = 0
  a_bar2  <- mean(curr_a)
  curr_a  <- curr_a - a_bar2
  curr_mu <- curr_mu + a_bar2

  # Save scalar traces
  samples_a[g, ]        <- curr_a
  samples_sigma2[g]     <- curr_sigma2
  samples_rho[g]        <- curr_rho
  samples_sigma_beta[g] <- curr_sigma_beta

  # Accumulate posterior field summaries after burn-in
  if (g > burnin) {
    mu_sum    <- mu_sum + curr_mu
    mu_sq_sum <- mu_sq_sum + curr_mu^2

    f_targ <- curr_a[target_year_idx] + as.vector(Phi %*% (curr_mu + curr_beta[, target_year_idx]))
    f_target_sum    <- f_target_sum + f_targ
    f_target_sq_sum <- f_target_sq_sum + f_targ^2
  }

  setTxtProgressBar(pb, g)
}

close(pb)
elapsed_time <- Sys.time() - start_time
cat("\nElapsed Time:\n")
print(elapsed_time)

# ------------------------------------------------------------------------------
# 6. POSTERIOR DIAGNOSTICS & SUMMARIES
# ------------------------------------------------------------------------------
idx_post <- (burnin + 1L):G

acc_rate <- mean(mh_accepted[idx_post])
rho_km_post <- samples_rho[idx_post] * scale_km
sig_eps_post <- sqrt(samples_sigma2[idx_post])
sig_beta_post <- samples_sigma_beta[idx_post]
kappa_post <- sqrt(8) / samples_rho[idx_post]

cat("\n=== POSTERIOR SUMMARIES (REAL DATA: ROME SQUARE) ===\n")
cat(sprintf("MH acceptance rate (post burn-in): %.1f%%\n", 100 * acc_rate))
cat(sprintf("sigma_eps  (degC): mean = %.4f | 95%% CI = [%.4f, %.4f]\n",
            mean(sig_eps_post), quantile(sig_eps_post, 0.025), quantile(sig_eps_post, 0.975)))
cat(sprintf("sigma_beta (degC): mean = %.4f | 95%% CI = [%.4f, %.4f]\n",
            mean(sig_beta_post), quantile(sig_beta_post, 0.025), quantile(sig_beta_post, 0.975)))
cat(sprintf("rho        (km)  : mean = %.1f   | 95%% CI = [%.1f, %.1f]\n",
            mean(rho_km_post), quantile(rho_km_post, 0.025), quantile(rho_km_post, 0.975)))
cat(sprintf("kappa      (unit): mean = %.4f | 95%% CI = [%.4f, %.4f]\n",
            mean(kappa_post), quantile(kappa_post, 0.025), quantile(kappa_post, 0.975)))

mu_post_mean <- mu_sum / n_post
f_targ_mean  <- f_target_sum / n_post
rmse_targ    <- sqrt(mean((Y_mat[, target_year_idx] - f_targ_mean)^2))
cat(sprintf("In-sample RMSE (%d): %.4f degC\n", years[target_year_idx], rmse_targ))

# --- Plot 1: Hyperparameter Traceplots ---
par(mfrow = c(2, 2))
plot(sqrt(samples_sigma2), type = "l", main = "sigma_eps (degC)", ylab = "degC")
abline(v = burnin, col = "red", lty = 2)

plot(samples_sigma_beta, type = "l", main = "sigma_beta (degC)", ylab = "degC")
abline(v = burnin, col = "red", lty = 2)

plot(samples_rho * scale_km, type = "l", main = "Spatial range rho (km)", ylab = "km")
abline(v = burnin, col = "red", lty = 2)

# Global yearly effect a_i (1983-2007) with 95% credible bands
a_mean <- colMeans(samples_a[idx_post, ])
a_lo   <- apply(samples_a[idx_post, ], 2, quantile, probs = 0.025)
a_hi   <- apply(samples_a[idx_post, ], 2, quantile, probs = 0.975)
plot(years, a_mean, type = "b", pch = 16, ylim = range(c(a_lo, a_hi)),
     main = "Global yearly effect a_i (1983-2007)", ylab = "degC", xlab = "Year")
arrows(years, a_lo, years, a_hi, angle = 90, code = 3, length = 0.04, col = "gray40")
abline(h = 0, lty = 2, col = "red")

# --- Plot 2: Spatial mean field mu(s) & 2007 reconstruction ---
par(mfrow = c(1, 2))
pal <- hcl.colors(100, "RdYlBu", rev = TRUE)

mu_obs <- as.vector(Phi %*% mu_post_mean)
z_rng  <- range(c(mu_obs, Y_mat[, target_year_idx], f_targ_mean))
col_idx <- function(v) pmax(1, pmin(100, round(1 + 99 * (v - z_rng[1]) / diff(z_rng))))

plot(real_data$coords, col = pal[col_idx(mu_obs)], pch = 15, cex = 1.4, asp = 1,
     main = "Posterior mean field mu(s) [degC]", xlab = "x (unit)", ylab = "y (unit)")
polygon(real_data$boundary, border = "black", lwd = 1.5)

plot(real_data$coords, col = pal[col_idx(f_targ_mean)], pch = 15, cex = 1.4, asp = 1,
     main = sprintf("Posterior mean %d [RMSE=%.2f C]", years[target_year_idx], rmse_targ),
     xlab = "x (unit)", ylab = "y (unit)")
polygon(real_data$boundary, border = "black", lwd = 1.5)
