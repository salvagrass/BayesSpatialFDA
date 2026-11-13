# SEQUENZA DI COMANDI DA LANCIARE IN CONSOLE (NON LANCIARE LO SCRIPT)
# rm(list = ls())
# setwd("~/Desktop/Progetto_Bayesian/BayesSpatialFDA/experiments")
# library(Matrix)
# library(MASS)
# library(Rcpp)
# library(RcppEigen)
# Rcpp::sourceCpp("../R/nngp_cpp.cpp")
# source("../R/nngp_functions.R")
# source("../R/nngp_marginal.R")
# source("test_nngp_marginal.R")

# ==============================================================================
# PROGETTO BAYESIANO: NNGP IBRIDO GIBBS/METROPOLIS (COMPLETO DI TUTTE LE TAVOLE)
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
Sigma_true <- sigma2_spaziale_true * exp(-phi_true * as.matrix(dist(coords)))
mu_true <- MASS::mvrnorm(1, mu = rep(15.0, N), Sigma = Sigma_true)

beta_true_raw <- matrix(NA, nrow = N, ncol = M)
Y_mat <- matrix(NA, nrow = N, ncol = M)
for (i in 1:M) {
  beta_true_raw[, i] <- MASS::mvrnorm(1, mu = rep(0, N), Sigma = Sigma_true)
  Y_mat[, i] <- mu_true + beta_true_raw[, i] + rnorm(N, 0, sqrt(sigma2_errore_true))
}

# Anomalia vera centrata per identificabilità rispetto alla media locale mu_i
beta_true <- beta_true_raw - rowMeans(beta_true_raw)

nngp_struct <- build_nngp_neighbors(coords, m = 10)

# SETUP MCMC
n_iter <- 50000
burn_in <- 2000
thin_step <- 20

samples_mu <- matrix(0, nrow = n_iter, ncol = N)
samples_sigma2_err <- numeric(n_iter)
samples_sigma2_spaz <- numeric(n_iter)
samples_phi <- numeric(n_iter)

# Accumulatori Online (Welford) per Media e Varianza a Posteriori delle Anomalie Beta
beta_post_mean <- matrix(0, nrow = N, ncol = M)
beta_M2        <- matrix(0, nrow = N, ncol = M)
n_thin_beta    <- 0

curr_mu <- rep(mean(Y_mat), N)
curr_beta <- matrix(0, nrow = N, ncol = M)
curr_sigma2_err <- 1.0
curr_sigma2_spaz <- 0.5
curr_phi <- 0.5
prior_mu_prec_scalar <- 1e-4

R_phi <- build_nngp_precision(nngp_struct, 1.0, curr_phi)
ldet_R_curr <- attr(R_phi, "log_det")

prop_var_phi <- 0.05
n_acc_phi <- 0
n_acc_total <- 0

cat("\nAvvio Sampler Ibrido Vettorizzato (", n_iter, " iterazioni)...\n")
pb <- txtProgressBar(min = 1, max = n_iter, style = 3)
start_time <- Sys.time()

for (iter in 1:n_iter) {

  # STEP 1: mu
  prec_post_mu <- (M / curr_sigma2_err) + prior_mu_prec_scalar
  b_post_mu <- rowSums(Y_mat - curr_beta) / curr_sigma2_err
  curr_mu <- rnorm(N, mean = b_post_mu / prec_post_mu, sd = sqrt(1 / prec_post_mu))

  # STEP 2: beta (Risoluzione C++ SimplicialLLT con fattorizzazione simbolica fissa)
  Z_norm <- matrix(rnorm(N * M), nrow = N, ncol = M)
  curr_beta <- sample_beta_cpp(R_phi, curr_sigma2_spaz, curr_sigma2_err, Y_mat - curr_mu, Z_norm)

  # STEP 3: Varianza Errore
  eps_mat <- Y_mat - curr_mu - curr_beta
  curr_sigma2_err <- 1 / rgamma(1, shape = (N * M) / 2 + 0.01, rate = sum(eps_mat^2) / 2 + 0.01)

  # STEP 4: Varianza Spaziale (Gibbs Coniugato)
  quad_R_phi <- quad_form_cpp(R_phi, curr_beta)
  shape_spaz <- (N * M) / 2 + 0.01
  rate_spaz  <- quad_R_phi / 2 + 0.01
  curr_sigma2_spaz <- 1 / rgamma(1, shape = shape_spaz, rate = rate_spaz)

  # STEP 5: METROPOLIS Lazy per PHI
  log_phi_prop <- rnorm(1, mean = log(curr_phi), sd = sqrt(prop_var_phi))
  phi_prop <- exp(log_phi_prop)

  cpp_prop <- eval_nngp_proposal_cpp(
    N = as.integer(N),
    coords_ord = nngp_struct$coords_ord,
    nn_idx = nngp_struct$nn_idx,
    ord = nngp_struct$ord_indices,
    sigma2 = 1.0,
    phi = phi_prop,
    beta = curr_beta
  )

  ldet_curr <- ldet_R_curr - N * log(curr_sigma2_spaz)
  ldet_prop <- cpp_prop$log_det - N * log(curr_sigma2_spaz)

  quad_curr <- quad_R_phi / curr_sigma2_spaz
  quad_prop <- cpp_prop$quad_form / curr_sigma2_spaz

  ll_curr <- 0.5 * M * ldet_curr - 0.5 * quad_curr
  ll_prop <- 0.5 * M * ldet_prop - 0.5 * quad_prop

  if (log(runif(1)) < (ll_prop + dnorm(log_phi_prop, 0, 5, log = TRUE) - ll_curr - dnorm(log(curr_phi), 0, 5, log = TRUE))) {
    curr_phi <- phi_prop
    R_phi <- assemble_nngp_precision(cpp_prop, N)
    ldet_R_curr <- cpp_prop$log_det
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

  # Aggiornamento Online di Media e Varianza dell'Anomalia Centrata
  if (iter > burn_in && (iter - burn_in) %% thin_step == 0) {
    n_thin_beta <- n_thin_beta + 1
    curr_beta_anom <- curr_beta - rowMeans(curr_beta)
    delta_beta <- curr_beta_anom - beta_post_mean
    beta_post_mean <- beta_post_mean + delta_beta / n_thin_beta
    beta_M2 <- beta_M2 + delta_beta * (curr_beta_anom - beta_post_mean)
  }

  setTxtProgressBar(pb, iter)
}
close(pb)
elapsed_sec <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
cat(sprintf("\nTempo totale: %.2f secondi (%.2f minuti)\n", elapsed_sec, elapsed_sec / 60))

# Matrici N x M di Varianza e Deviazione Standard a posteriori delle anomalie Beta
beta_post_var <- (beta_M2 / (n_thin_beta - 1)) * (M / (M - 1))
beta_post_sd  <- sqrt(beta_post_var)
beta_node_sd  <- sqrt(rowMeans(beta_post_var))

# ==============================================================================
# DIAGNOSTICA NUMERICHE POST-MCMC
# ==============================================================================
idx_thin <- seq(burn_in + 1, n_iter, by = thin_step)
samples_sigma2_err_thin <- samples_sigma2_err[idx_thin]
samples_sigma2_spaz_thin <- samples_sigma2_spaz[idx_thin]
samples_phi_thin <- samples_phi[idx_thin]
Y_centered_obs <- Y_mat - rowMeans(Y_mat)

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

rmse_beta  <- sqrt(mean((beta_post_mean - beta_true)^2))
cor_beta   <- cor(as.vector(beta_post_mean), as.vector(beta_true))
ci95_low   <- beta_post_mean - 1.96 * beta_post_sd
ci95_high  <- beta_post_mean + 1.96 * beta_post_sd
cov95_beta <- 100 * mean(beta_true >= ci95_low & beta_true <= ci95_high)

cat("\n====================================================================\n")
cat(" RIASSUNTO POSTERIORE E DIAGNOSTICA MCMC (N_thin =", length(idx_thin), ")\n")
cat("====================================================================\n")
print(summary_tab, row.names = FALSE, digits = 4)
cat(sprintf("Tasso di accettazione post burn-in (phi): %.2f%%\n", 100 * n_acc_total / (n_iter - burn_in)))
cat("--------------------------------------------------------------------\n")
cat(" ANALISI A POSTERIORI DEL CAMPO LATENTE BETA (Matrice N x M = 400 x 25)\n")
cat("--------------------------------------------------------------------\n")
cat(sprintf("Correlazione col Vero  : R = %.4f | RMSE Globale = %.4f\n", cor_beta, rmse_beta))
cat(sprintf("Varianza Var(Beta | Y) : Media = %.4f | Min (centro) = %.4f | Max (bordi) = %.4f\n",
            mean(beta_post_var), min(beta_node_sd^2), max(beta_node_sd^2)))
cat(sprintf("Dev. Std SD(Beta | Y)  : Media = %.4f | Min (centro) = %.4f | Max (bordi) = %.4f\n",
            mean(beta_node_sd), min(beta_node_sd), max(beta_node_sd)))
cat(sprintf("Copertura IC 95%%       : %.2f%%\n", cov95_beta))
cat("====================================================================\n")

# ==============================================================================
# TAVOLA GRAFICA 1 (3x2): TRACEPLOT E DENSITA' PARAMETRI GLOBALI
# ==============================================================================
par(mfrow = c(3, 2), pty = "m", mar = c(4, 4, 3, 1))

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

# ==============================================================================
# TAVOLA GRAFICA 2 (1x2): CURVA TEMPORALE BETA E RECUPERO GLOBALE
# ==============================================================================
par(mfrow = c(1, 2), pty = "m", mar = c(4, 4, 3, 1.5))

sito_centro <- which.min((coords[, 1] - 2.5)^2 + (coords[, 2] - 2.5)^2)
y_centr <- Y_centered_obs[sito_centro, ]
m_c <- beta_post_mean[sito_centro, ]
s_c <- beta_post_sd[sito_centro, ]
ylim_c <- range(c(y_centr, beta_true[sito_centro, ], m_c - 2 * s_c, m_c + 2 * s_c))

plot(1:M, beta_true[sito_centro, ], type = "l", lwd = 2, col = "black", ylim = ylim_c,
     xlab = "Anno / Tempo (j = 1..M)", ylab = expression(beta[ij]), cex.main = 0.95,
     main = paste0("Curva Beta - Nodo Centrale (SD = ", round(mean(s_c), 3), ")"))
polygon(c(1:M, M:1), c(m_c + 2 * s_c, rev(m_c - 2 * s_c)), col = rgb(0.1, 0.5, 0.8, 0.20), border = NA)
polygon(c(1:M, M:1), c(m_c + s_c, rev(m_c - s_c)), col = rgb(0.1, 0.5, 0.8, 0.35), border = NA)
lines(1:M, m_c, col = "dodgerblue4", lwd = 2, lty = 2)
points(1:M, y_centr, pch = 16, cex = 0.7, col = "firebrick")
legend("topright", legend = c("Beta Vero", "E[Beta|Y]", "+/- 1 SD", "+/- 2 SD", "Y centrato"),
       col = c("black", "dodgerblue4", rgb(0.1, 0.5, 0.8, 0.45), rgb(0.1, 0.5, 0.8, 0.25), "firebrick"),
       lty = c(1, 2, NA, NA, NA), pch = c(NA, NA, 15, 15, 16), lwd = 2, bty = "n", cex = 0.75)

plot(as.vector(beta_true), as.vector(beta_post_mean), pch = 16, cex = 0.3, col = rgb(0, 0.4, 0.8, 0.3),
     xlab = "True Anomaly (Beta Vero)", ylab = "Posterior Mean E[Beta | Y]", cex.main = 0.95,
     main = paste0("Recupero Globale Anomalie Beta (R = ", round(cor_beta, 3), ")"))
abline(0, 1, col = "red", lwd = 2)

# ==============================================================================
# TAVOLA GRAFICA 3 (2x2): GRIGLIA 20x20 AUTENTICA SENZA BORDI BIANCHI
# ==============================================================================
x_grid <- seq(0, 5, length.out = n_side)
y_grid <- seq(0, 5, length.out = n_side)
pal_spaz <- hcl.colors(100, "Spectral", rev = TRUE)

par(mfrow = c(2, 2), pty = "s", mar = c(3.8, 3.8, 2.8, 1.2))

zlim_anom <- range(c(beta_true[, 1], beta_post_mean[, 1]))

# 1. Anomalia Vera (t = 1)
image(x_grid, y_grid, matrix(beta_true[, 1], nrow = n_side, ncol = n_side),
      col = pal_spaz, zlim = zlim_anom, xlim = c(0, 5), ylim = c(0, 5),
      xaxs = "i", yaxs = "i", xlab = "x", ylab = "y", main = "True anomaly (t = 1)")
box()

# 2. Media a Posteriori NNGP (t = 1)
image(x_grid, y_grid, matrix(beta_post_mean[, 1], nrow = n_side, ncol = n_side),
      col = pal_spaz, zlim = zlim_anom, xlim = c(0, 5), ylim = c(0, 5),
      xaxs = "i", yaxs = "i", xlab = "x", ylab = "y", main = "Posterior mean NNGP (t = 1)")
box()

# 3. Deviazione Standard a Posteriori NNGP
image(x_grid, y_grid, matrix(beta_node_sd, nrow = n_side, ncol = n_side),
      col = pal_spaz, xlim = c(0, 5), ylim = c(0, 5),
      xaxs = "i", yaxs = "i", xlab = "x", ylab = "y", main = "Posterior SD NNGP")
box()

# 4. Errore di Stima (Media a Posteriori - Vero, t = 1)
err_t1 <- beta_post_mean[, 1] - beta_true[, 1]
image(x_grid, y_grid, matrix(err_t1, nrow = n_side, ncol = n_side),
      col = pal_spaz, xlim = c(0, 5), ylim = c(0, 5),
      xaxs = "i", yaxs = "i", xlab = "x", ylab = "y", main = "Error (Posterior mean - Truth)")
box()

par(mfrow = c(1, 1), pty = "m")
