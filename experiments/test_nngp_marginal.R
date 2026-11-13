# ==============================================================================
# PROGETTO BAYESIANO: NNGP IBRIDO GIBBS/METROPOLIS (ALTA VELOCITA')
# ==============================================================================

library(Matrix)
library(MASS)

# SETUP DATI
set.seed(42)
cat("\nGenerazione dati su griglia 20x20...\n")
n_side <- 20
coords <- as.matrix(expand.grid(x = seq(0, 5, length.out = n_side),
                                y = seq(0, 5, length.out = n_side)))
N <- nrow(coords)
M <- 25

sigma2_spaziale_true <- 2.0
phi_true <- 1.2
sigma2_errore_true <- 0.1
mu_true <- MASS::mvrnorm(1, mu = rep(15.0, N), Sigma = sigma2_spaziale_true * exp(-phi_true * as.matrix(dist(coords))))

Y_mat <- matrix(NA, nrow = N, ncol = M)
for (i in 1:M) {
  Y_mat[, i] <- mu_true + MASS::mvrnorm(1, mu = rep(0, N), Sigma = sigma2_spaziale_true * exp(-phi_true * as.matrix(dist(coords)))) + rnorm(N, 0, sqrt(sigma2_errore_true))
}

nngp_struct <- build_nngp_neighbors(coords, m = 10)

# SETUP MCMC
n_iter <- 50000
burn_in <- 2000
thin_step <- 20

samples_mu <- matrix(0, nrow = n_iter, ncol = N)
samples_sigma2_err <- numeric(n_iter)
samples_sigma2_spaz <- numeric(n_iter)
samples_phi <- numeric(n_iter)

curr_mu <- rep(mean(Y_mat), N)
curr_beta <- matrix(0, nrow = N, ncol = M)
curr_sigma2_err <- 1.0
curr_sigma2_spaz <- 0.5
curr_phi <- 0.5
prior_mu_prec <- Diagonal(N, 1e-4)

R_phi <- build_nngp_precision(nngp_struct, 1.0, curr_phi)
ldet_R_curr <- attr(R_phi, "log_det")
Q_nngp <- R_phi / curr_sigma2_spaz

prop_var_phi <- 0.05
n_acc_phi <- 0
n_acc_total <- 0

cat("\nAvvio Sampler Ibrido Vettorizzato (", n_iter, " iterazioni)...\n")
pb <- txtProgressBar(min = 1, max = n_iter, style = 3)
start_time <- Sys.time()

for (iter in 1:n_iter) {

  # STEP 1: mu
  Q_post_mu <- Diagonal(N, M / curr_sigma2_err) + prior_mu_prec
  sum_residui_mu <- rowSums(Y_mat - curr_beta)
  b_post_mu <- sum_residui_mu / curr_sigma2_err
  curr_mu <- rnorm(N, mean = as.vector(b_post_mu / diag(Q_post_mu)), sd = sqrt(1 / diag(Q_post_mu)))

  # STEP 2: beta (Vettorizzato su dsCMatrix simmetrica)
  Q_post_beta <- Q_nngp + Diagonal(N, 1 / curr_sigma2_err)
  chol_Q_beta <- chol(Q_post_beta)
  b_post_beta_mat <- (Y_mat - curr_mu) / curr_sigma2_err
  mean_beta_post_mat <- as.matrix(solve(chol_Q_beta, solve(t(chol_Q_beta), b_post_beta_mat)))
  curr_beta <- as.matrix(mean_beta_post_mat + solve(chol_Q_beta, matrix(rnorm(N * M), nrow = N, ncol = M)))

  # STEP 3: Varianza Errore
  eps_mat <- Y_mat - curr_mu - curr_beta
  curr_sigma2_err <- 1 / rgamma(1, shape = (N * M) / 2 + 0.01, rate = sum(eps_mat^2) / 2 + 0.01)

  # STEP 4: Varianza Spaziale (Gibbs Coniugato)
  quad_R_phi <- sum(curr_beta * (R_phi %*% curr_beta))
  shape_spaz <- (N * M) / 2 + 0.01
  rate_spaz  <- quad_R_phi / 2 + 0.01
  curr_sigma2_spaz <- 1 / rgamma(1, shape = shape_spaz, rate = rate_spaz)
  Q_nngp <- R_phi / curr_sigma2_spaz

  # STEP 5: METROPOLIS per PHI (Con log-determinante analitico O(N) gratuito!)
  log_phi_prop <- rnorm(1, mean = log(curr_phi), sd = sqrt(prop_var_phi))
  phi_prop <- exp(log_phi_prop)

  R_phi_prop <- build_nngp_precision(nngp_struct, 1.0, phi_prop)
  ldet_R_prop <- attr(R_phi_prop, "log_det")

  # Nota: il termine -N*log(curr_sigma2_spaz) si cancella esattamente tra curr e prop,
  # ma lo includiamo per identità numerica assoluta col tuo script originale
  ldet_curr <- ldet_R_curr - N * log(curr_sigma2_spaz)
  ldet_prop <- ldet_R_prop - N * log(curr_sigma2_spaz)

  quad_curr <- quad_R_phi / curr_sigma2_spaz
  quad_prop <- sum(curr_beta * (R_phi_prop %*% curr_beta)) / curr_sigma2_spaz

  ll_curr <- 0.5 * M * ldet_curr - 0.5 * quad_curr
  ll_prop <- 0.5 * M * ldet_prop - 0.5 * quad_prop

  if (log(runif(1)) < (ll_prop + dnorm(log_phi_prop, 0, 5, log = TRUE) - ll_curr - dnorm(log(curr_phi), 0, 5, log = TRUE))) {
    curr_phi <- phi_prop
    R_phi <- R_phi_prop
    ldet_R_curr <- ldet_R_prop
    Q_nngp <- R_phi / curr_sigma2_spaz
    n_acc_phi <- n_acc_phi + 1
    if (iter > burn_in) n_acc_total <- n_acc_total + 1
  }

  if (iter <= burn_in && iter %% 100 == 0) {
    if (n_acc_phi / 100 > 0.44) prop_var_phi <- prop_var_phi * 1.2
    if (n_acc_phi / 100 < 0.23) prop_var_phi <- prop_var_phi * 0.8
    n_acc_phi <- 0
  }

  samples_mu[iter, ] <- curr_mu
  samples_sigma2_err[iter] <- curr_sigma2_err
  samples_sigma2_spaz[iter] <- curr_sigma2_spaz
  samples_phi[iter] <- curr_phi

  setTxtProgressBar(pb, iter)
}
close(pb)
elapsed_sec <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
cat(sprintf("\nTempo totale: %.2f secondi (%.2f minuti)\n", elapsed_sec, elapsed_sec / 60))

# ==============================================================================
# DIAGNOSTICA NUMERICHE E GRAFICHE POST-MCMC
# ==============================================================================
idx_thin <- seq(burn_in + 1, n_iter, by = thin_step)
samples_sigma2_err_thin <- samples_sigma2_err[idx_thin]
samples_sigma2_spaz_thin <- samples_sigma2_spaz[idx_thin]
samples_phi_thin <- samples_phi[idx_thin]

calc_ess <- function(x) {
  n <- length(x)
  ac <- acf(x, lag.max = min(100, n - 1), plot = FALSE)$acf[-1]
  pos <- which(ac < 0.05)
  k <- if (length(pos) > 0) pos[1] - 1 else length(ac)
  if (k == 0) return(n)
  return(n / (1 + 2 * sum(ac[1:k])))
}

summary_tab <- data.frame(
  Parametro = c("sigma2_spaz", "phi", "sigma2_err"),
  Vero      = c(sigma2_spaziale_true, phi_true, sigma2_errore_true),
  Media     = c(mean(samples_sigma2_spaz_thin), mean(samples_phi_thin), mean(samples_sigma2_err_thin)),
  Mediana   = c(median(samples_sigma2_spaz_thin), median(samples_phi_thin), median(samples_sigma2_err_thin)),
  CI_2.5    = c(quantile(samples_sigma2_spaz_thin, 0.025), quantile(samples_phi_thin, 0.025), quantile(samples_sigma2_err_thin, 0.025)),
  CI_97.5   = c(quantile(samples_sigma2_spaz_thin, 0.975), quantile(samples_phi_thin, 0.975), quantile(samples_sigma2_err_thin, 0.975)),
  ESS       = round(c(calc_ess(samples_sigma2_spaz_thin), calc_ess(samples_phi_thin), calc_ess(samples_sigma2_err_thin)), 1)
)

cat("\n====================================================================\n")
cat(" RIASSUNTO POSTERIORE E DIAGNOSTICA MCMC (N_thin =", length(idx_thin), ")\n")
cat("====================================================================\n")
print(summary_tab, row.names = FALSE, digits = 4)
cat(sprintf("Tasso di accettazione post burn-in (phi): %.2f%%\n", 100 * n_acc_total / (n_iter - burn_in)))
cat("====================================================================\n")

# PLOT E DIAGNOSTICA
par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))

ylim_spaz <- range(c(samples_sigma2_spaz_thin, sigma2_spaziale_true))
plot(idx_thin, samples_sigma2_spaz_thin, type = "l", col = "darkorange",
     xlab = "Iterazione", ylab = "sigma2_spaz", main = "Traceplot Spaziale", ylim = ylim_spaz)
abline(h = sigma2_spaziale_true, col = "red", lwd = 2)

plot(density(samples_sigma2_spaz_thin), col = "darkorange", lwd = 2, main = "Densità a Posteriori (sigma2_spaz)")
abline(v = sigma2_spaziale_true, col = "red", lwd = 2, lty = 2)
abline(v = quantile(samples_sigma2_spaz_thin, probs = c(0.025, 0.975)), col = "darkorange", lty = 3)

ylim_phi <- range(c(samples_phi_thin, phi_true))
plot(idx_thin, samples_phi_thin, type = "l", col = "cyan4",
     xlab = "Iterazione", ylab = "phi", main = "Traceplot Decadimento", ylim = ylim_phi)
abline(h = phi_true, col = "red", lwd = 2)

plot(density(samples_phi_thin), col = "cyan4", lwd = 2, main = "Densità a Posteriori (phi)")
abline(v = phi_true, col = "red", lwd = 2, lty = 2)
abline(v = quantile(samples_phi_thin, probs = c(0.025, 0.975)), col = "cyan4", lty = 3)

ylim_err <- range(c(samples_sigma2_err_thin, sigma2_errore_true))
plot(idx_thin, samples_sigma2_err_thin, type = "l", col = "grey30",
     xlab = "Iterazione", ylab = "sigma2_err", main = "Traceplot Errore", ylim = ylim_err)
abline(h = sigma2_errore_true, col = "red", lwd = 2)

plot(density(samples_sigma2_err_thin), col = "grey30", lwd = 2, main = "Densità a Posteriori (sigma2_err)")
abline(v = sigma2_errore_true, col = "red", lwd = 2, lty = 2)
abline(v = quantile(samples_sigma2_err_thin, probs = c(0.025, 0.975)), col = "grey30", lty = 3)

par(mfrow = c(1, 1))
