# ==============================================================================
# PROGETTO BAYESIANO: NNGP LATENTE VETTORIZZATO (ULTRA-VELOCE)
# ==============================================================================

library(Matrix)
library(MASS)
library(devtools)

load_all("R")

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

# SETUP MCMC: Grazie alla vettorizzazione faremo 25000 iterazioni in pochi secondi!
n_iter <- 25000
burn_in <- 5000
thin_step <- 40

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

# Varianze Proposte Iniziali per Metropolis
prop_var_sig2 <- 0.05
prop_var_phi <- 0.05
n_acc_sig2 <- 0; n_acc_phi <- 0

Q_nngp <- build_nngp_precision(nngp_struct, curr_sigma2_spaz, curr_phi)

cat("\nAvvio Sampler VETTORIZZATO (", n_iter, "iterazioni )...\n")
pb <- txtProgressBar(min = 1, max = n_iter, style = 3)
start_time <- Sys.time()

for (iter in 1:n_iter) {

  # -----------------------------------------------------------------
  # STEP 1: Aggiornamento mu (Vettorizzato)
  # -----------------------------------------------------------------
  Q_post_mu <- Diagonal(N, M / curr_sigma2_err) + prior_mu_prec
  # Sostituiamo il ciclo for con rowSums (Somma istantanea per riga)
  sum_residui_mu <- rowSums(Y_mat - curr_beta)
  b_post_mu <- sum_residui_mu / curr_sigma2_err
  curr_mu <- rnorm(N, mean = as.vector(b_post_mu / diag(Q_post_mu)), sd = sqrt(1 / diag(Q_post_mu)))

  # -----------------------------------------------------------------
  # STEP 2: Aggiornamento beta_tilde (Vettorizzato con matrici M)
  # -----------------------------------------------------------------
  Q_post_beta <- Q_nngp + Diagonal(N, 1 / curr_sigma2_err)
  chol_Q_beta <- chol(Q_post_beta)

  # Invece di risolvere 25 sistemi uno alla volta, R risolve la matrice intera N x M in un colpo solo!
  b_post_beta_mat <- (Y_mat - curr_mu) / curr_sigma2_err
  mean_beta_post_mat <- as.matrix(solve(chol_Q_beta, solve(t(chol_Q_beta), b_post_beta_mat)))

  # Estrazione del rumore per tutta la matrice
  Z_noise <- matrix(rnorm(N * M), nrow = N, ncol = M)
  curr_beta <- as.matrix(mean_beta_post_mat + solve(chol_Q_beta, Z_noise))

  # -----------------------------------------------------------------
  # STEP 3: Varianza Errore (Vettorizzato)
  # -----------------------------------------------------------------
  eps_mat <- Y_mat - curr_mu - curr_beta
  curr_sigma2_err <- 1 / rgamma(1, shape = (N * M) / 2 + 0.01, rate = sum(eps_mat^2) / 2 + 0.01)

  # -----------------------------------------------------------------
  # STEP 4: METROPOLIS-HASTINGS LATENTE (Forme quadratiche istantanee)
  # -----------------------------------------------------------------
  ldet_curr <- as.numeric(determinant(Q_nngp, logarithm = TRUE)$modulus)
  # Calcolo esatto traccia: sum(B * (Q %*% B)). Zero cicli for!
  quad_curr <- sum(curr_beta * (Q_nngp %*% curr_beta))
  ll_curr <- 0.5 * M * ldet_curr - 0.5 * quad_curr

  # Proposta Varianza Spaziale
  log_sig2_prop <- rnorm(1, mean = log(curr_sigma2_spaz), sd = sqrt(prop_var_sig2))
  sig2_prop <- exp(log_sig2_prop)
  Q_prop_sig <- build_nngp_precision(nngp_struct, sig2_prop, curr_phi)
  ldet_prop_sig <- as.numeric(determinant(Q_prop_sig, logarithm = TRUE)$modulus)
  quad_prop_sig <- sum(curr_beta * (Q_prop_sig %*% curr_beta))
  ll_prop_sig <- 0.5 * M * ldet_prop_sig - 0.5 * quad_prop_sig

  if (log(runif(1)) < (ll_prop_sig + dnorm(log_sig2_prop,0,5,log=T) - ll_curr - dnorm(log(curr_sigma2_spaz),0,5,log=T))) {
    curr_sigma2_spaz <- sig2_prop
    Q_nngp <- Q_prop_sig
    ll_curr <- ll_prop_sig
    n_acc_sig2 <- n_acc_sig2 + 1
  }

  # Proposta Decadimento (phi)
  log_phi_prop <- rnorm(1, mean = log(curr_phi), sd = sqrt(prop_var_phi))
  phi_prop <- exp(log_phi_prop)
  Q_prop_phi <- build_nngp_precision(nngp_struct, curr_sigma2_spaz, phi_prop)
  ldet_prop_phi <- as.numeric(determinant(Q_prop_phi, logarithm = TRUE)$modulus)
  quad_prop_phi <- sum(curr_beta * (Q_prop_phi %*% curr_beta))
  ll_prop_phi <- 0.5 * M * ldet_prop_phi - 0.5 * quad_prop_phi

  if (log(runif(1)) < (ll_prop_phi + dnorm(log_phi_prop,0,5,log=T) - ll_curr - dnorm(log(curr_phi),0,5,log=T))) {
    curr_phi <- phi_prop
    Q_nngp <- Q_prop_phi
    n_acc_phi <- n_acc_phi + 1
  }

  # Tuning Adattivo
  if (iter <= burn_in && iter %% 100 == 0) {
    if (n_acc_sig2/100 > 0.44) prop_var_sig2 <- prop_var_sig2 * 1.2
    if (n_acc_sig2/100 < 0.23) prop_var_sig2 <- prop_var_sig2 * 0.8
    if (n_acc_phi/100 > 0.44) prop_var_phi <- prop_var_phi * 1.2
    if (n_acc_phi/100 < 0.23) prop_var_phi <- prop_var_phi * 0.8
    n_acc_sig2 <- 0; n_acc_phi <- 0
  }

  samples_mu[iter, ] <- curr_mu
  samples_sigma2_err[iter] <- curr_sigma2_err
  samples_sigma2_spaz[iter] <- curr_sigma2_spaz
  samples_phi[iter] <- curr_phi

  setTxtProgressBar(pb, iter)
}
close(pb)
cat("\nTempo totale:", round(Sys.time() - start_time, 2), "secondi\n")

# ==============================================================================
# ESTRAZIONE CAMPIONI E DIAGNOSTICA VISIVA AVANZATA (ZOOM E DENSITA')
# ==============================================================================
idx_thin <- seq(burn_in + 1, n_iter, by = thin_step)

samples_sigma2_err_thin <- samples_sigma2_err[idx_thin]
samples_sigma2_spaz_thin <- samples_sigma2_spaz[idx_thin]
samples_phi_thin <- samples_phi[idx_thin]

# Impostiamo una griglia 3x2: Traceplot a sinistra, Densità a destra
par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))

# --- 1. VARIANZA SPAZIALE (sigma2_spaz) ---
# Calcolo automatico dello zoom (inclusivo del valore vero per non farlo sparire)
ylim_spaz <- range(c(samples_sigma2_spaz_thin, sigma2_spaziale_true))

plot(idx_thin, samples_sigma2_spaz_thin, type = "l", col = "darkorange",
     xlab = "Iterazione", ylab = "sigma2_spaz",
     main = "Traceplot Spaziale (Zoom Dinamico)", ylim = ylim_spaz)
abline(h = sigma2_spaziale_true, col = "red", lwd = 2)

plot(density(samples_sigma2_spaz_thin), col = "darkorange", lwd = 2,
     main = "Densità a Posteriori (sigma2_spaz)", xlab = "Valore", ylab = "Densità")
abline(v = sigma2_spaziale_true, col = "red", lwd = 2, lty = 2)
# Aggiungiamo l'intervallo di credibilità al 95%
abline(v = quantile(samples_sigma2_spaz_thin, probs = c(0.025, 0.975)), col = "darkorange", lty = 3)

# --- 2. DECADIMENTO SPAZIALE (phi) ---
ylim_phi <- range(c(samples_phi_thin, phi_true))

plot(idx_thin, samples_phi_thin, type = "l", col = "cyan4",
     xlab = "Iterazione", ylab = "phi",
     main = "Traceplot Decadimento (Zoom Dinamico)", ylim = ylim_phi)
abline(h = phi_true, col = "red", lwd = 2)

plot(density(samples_phi_thin), col = "cyan4", lwd = 2,
     main = "Densità a Posteriori (phi)", xlab = "Valore", ylab = "Densità")
abline(v = phi_true, col = "red", lwd = 2, lty = 2)
abline(v = quantile(samples_phi_thin, probs = c(0.025, 0.975)), col = "cyan4", lty = 3)

# --- 3. VARIANZA ERRORE (sigma2_err) ---
ylim_err <- range(c(samples_sigma2_err_thin, sigma2_errore_true))

plot(idx_thin, samples_sigma2_err_thin, type = "l", col = "grey30",
     xlab = "Iterazione", ylab = "sigma2_err",
     main = "Traceplot Errore (Zoom Dinamico)", ylim = ylim_err)
abline(h = sigma2_errore_true, col = "red", lwd = 2)

plot(density(samples_sigma2_err_thin), col = "grey30", lwd = 2,
     main = "Densità a Posteriori (sigma2_err)", xlab = "Valore", ylab = "Densità")
abline(v = sigma2_errore_true, col = "red", lwd = 2, lty = 2)
abline(v = quantile(samples_sigma2_err_thin, probs = c(0.025, 0.975)), col = "grey30", lty = 3)

par(mfrow = c(1, 1))
