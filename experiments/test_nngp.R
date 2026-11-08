# ==============================================================================
# PROGETTO BAYESIANO: GIBBS SAMPLER NNGP (N = 400, TUNING ADATTIVO)
# Modello: Y_i = mu + beta_tilde_i + epsilon_i
# ==============================================================================

library(Matrix)
library(MASS)
library(devtools)

# Assicurati di avere la funzione build_nngp_neighbors e build_nngp_precision in "R"
load_all("R")

# --- 1. SETUP E GENERAZIONE DEI DATI "GROUND TRUTH" ---
set.seed(42)
cat("\nGenerazione dati su griglia 20x20...\n")
n_side <- 20
coords <- as.matrix(expand.grid(x = seq(0, 5, length.out = n_side),
                                y = seq(0, 5, length.out = n_side)))
N <- nrow(coords)
M <- 25 # Anni di osservazione

# Parametri Veri
sigma2_spaziale_true <- 2.0
phi_true <- 1.2
sigma2_errore_true <- 0.1
media_globale_assoluta <- 15.0

# Generiamo la struttura spaziale
dist_mat <- as.matrix(dist(coords))
Sigma_true <- sigma2_spaziale_true * exp(-phi_true * dist_mat)

# Campo spaziale medio (Assorbe la temperatura assoluta)
mu_true <- MASS::mvrnorm(1, mu = rep(media_globale_assoluta, N), Sigma = Sigma_true)

# Generiamo i dati Y_i
Y_mat <- matrix(NA, nrow = N, ncol = M)
for (i in 1:M) {
  beta_i <- MASS::mvrnorm(1, mu = rep(0, N), Sigma = Sigma_true)
  Y_mat[, i] <- mu_true + beta_i + rnorm(N, mean = 0, sd = sqrt(sigma2_errore_true))
}

# --- 2. PRE-COMPUTAZIONE NNGP ---
nngp_struct <- build_nngp_neighbors(coords, m = 10)

# --- 3. SETUP MCMC ---
n_iter <- 5000
burn_in <- 1000

samples_mu <- matrix(0, nrow = n_iter, ncol = N)
samples_sigma2_err <- numeric(n_iter)
samples_sigma2_spaz <- numeric(n_iter)
samples_phi <- numeric(n_iter)

# Inizializzazione
curr_mu <- rep(mean(Y_mat), N)
curr_beta <- matrix(0, nrow = N, ncol = M)
curr_sigma2_err <- 1.0
curr_sigma2_spaz <- 0.5
curr_phi <- 0.5

Q_nngp <- build_nngp_precision(nngp_struct, curr_sigma2_spaz, curr_phi)
prior_mu_prec <- Diagonal(N, 1e-4)

# Variabili per il Tuning Adattivo del Metropolis-Hastings
prop_var_sig2 <- 0.05
prop_var_phi <- 0.05
n_acc_sig2 <- 0
n_acc_phi <- 0

cat("\nAvvio Gibbs Sampler NNGP con Tuning Adattivo su", N, "vertici...\n")
pb <- txtProgressBar(min = 1, max = n_iter, style = 3)
start_time <- Sys.time()

# --- 4. GIBBS SAMPLER ---
for (iter in 1:n_iter) {

  # ---------------------------------------------------------
  # STEP 1: Aggiornamento Campo Spaziale Globale (mu)
  # ---------------------------------------------------------
  Q_data <- Diagonal(N, M / curr_sigma2_err)
  Q_post_mu <- Q_data + prior_mu_prec

  sum_residui_mu <- rep(0, N)
  for (i in 1:M) {
    sum_residui_mu <- sum_residui_mu + (Y_mat[, i] - curr_beta[, i])
  }
  b_post_mu <- as.vector(sum_residui_mu / curr_sigma2_err)

  mean_mu_post <- as.vector(b_post_mu / diag(Q_post_mu))
  var_mu_post <- 1 / diag(Q_post_mu)
  curr_mu <- rnorm(N, mean = mean_mu_post, sd = sqrt(var_mu_post))

  # ---------------------------------------------------------
  # STEP 2: Aggiornamento Anomalie Spaziali Annuali (beta_tilde_i)
  # ---------------------------------------------------------
  Q_post_beta <- Q_nngp + Diagonal(N, 1 / curr_sigma2_err)
  chol_Q_beta <- chol(Q_post_beta)

  for (i in 1:M) {
    b_post_beta <- (Y_mat[, i] - curr_mu) / curr_sigma2_err
    mean_beta_post <- as.vector(solve(chol_Q_beta, solve(t(chol_Q_beta), b_post_beta)))
    z <- rnorm(N)
    curr_beta[, i] <- as.vector(mean_beta_post + solve(chol_Q_beta, z))
  }

  # ---------------------------------------------------------
  # STEP 3: Aggiornamento Varianza Errore (sigma2_err)
  # ---------------------------------------------------------
  eps_mat <- Y_mat - curr_mu - curr_beta
  shape_err <- (N * M) / 2 + 0.01
  rate_err <- sum(eps_mat^2) / 2 + 0.01
  curr_sigma2_err <- 1 / rgamma(1, shape = shape_err, rate = rate_err)

  # ---------------------------------------------------------
  # STEP 4: Metropolis-Hastings Sequenziale con Tuning
  # ---------------------------------------------------------

  # Calcoliamo la verosimiglianza dello stato attuale una sola volta
  ldet_curr <- as.numeric(determinant(Q_nngp, logarithm = TRUE)$modulus)
  quad_curr <- sum(sapply(1:M, function(i) as.numeric(t(curr_beta[,i]) %*% Q_nngp %*% curr_beta[,i])))
  ll_curr <- 0.5 * M * ldet_curr - 0.5 * quad_curr

  # A) Proposta e Aggiornamento per Varianza Spaziale (sigma2_spaz)
  log_sig2_prop <- rnorm(1, mean = log(curr_sigma2_spaz), sd = sqrt(prop_var_sig2))
  sig2_prop <- exp(log_sig2_prop)

  Q_prop_sig <- build_nngp_precision(nngp_struct, sig2_prop, curr_phi)
  ldet_prop_sig <- as.numeric(determinant(Q_prop_sig, logarithm = TRUE)$modulus)
  quad_prop_sig <- sum(sapply(1:M, function(i) as.numeric(t(curr_beta[,i]) %*% Q_prop_sig %*% curr_beta[,i])))
  ll_prop_sig <- 0.5 * M * ldet_prop_sig - 0.5 * quad_prop_sig

  lprior_curr_sig <- dnorm(log(curr_sigma2_spaz), 0, 5, log = TRUE)
  lprior_prop_sig <- dnorm(log_sig2_prop, 0, 5, log = TRUE)

  if (log(runif(1)) < (ll_prop_sig + lprior_prop_sig - ll_curr - lprior_curr_sig)) {
    curr_sigma2_spaz <- sig2_prop
    Q_nngp <- Q_prop_sig
    ll_curr <- ll_prop_sig # Aggiorniamo la likelihood base per il prossimo step
    n_acc_sig2 <- n_acc_sig2 + 1
  }

  # B) Proposta e Aggiornamento per Decadimento Spaziale (phi)
  log_phi_prop <- rnorm(1, mean = log(curr_phi), sd = sqrt(prop_var_phi))
  phi_prop <- exp(log_phi_prop)

  Q_prop_phi <- build_nngp_precision(nngp_struct, curr_sigma2_spaz, phi_prop)
  ldet_prop_phi <- as.numeric(determinant(Q_prop_phi, logarithm = TRUE)$modulus)
  quad_prop_phi <- sum(sapply(1:M, function(i) as.numeric(t(curr_beta[,i]) %*% Q_prop_phi %*% curr_beta[,i])))
  ll_prop_phi <- 0.5 * M * ldet_prop_phi - 0.5 * quad_prop_phi

  lprior_curr_phi <- dnorm(log(curr_phi), 0, 5, log = TRUE)
  lprior_prop_phi <- dnorm(log_phi_prop, 0, 5, log = TRUE)

  if (log(runif(1)) < (ll_prop_phi + lprior_prop_phi - ll_curr - lprior_curr_phi)) {
    curr_phi <- phi_prop
    Q_nngp <- Q_prop_phi
    n_acc_phi <- n_acc_phi + 1
  }

  # C) Tuning Adattivo (Ogni 50 iterazioni, SOLO durante il burn-in)
  if (iter <= burn_in && iter %% 50 == 0) {
    acc_rate_sig2 <- n_acc_sig2 / 50
    acc_rate_phi <- n_acc_phi / 50

    # Aggiusta la varianza della proposta per mantenere l'accettazione tra 23% e 44%
    if (acc_rate_sig2 > 0.44) prop_var_sig2 <- prop_var_sig2 * 1.5
    if (acc_rate_sig2 < 0.23) prop_var_sig2 <- prop_var_sig2 * 0.5

    if (acc_rate_phi > 0.44) prop_var_phi <- prop_var_phi * 1.5
    if (acc_rate_phi < 0.23) prop_var_phi <- prop_var_phi * 0.5

    # Reset contatori
    n_acc_sig2 <- 0
    n_acc_phi <- 0
  }

  # Salvataggio
  samples_mu[iter, ] <- curr_mu
  samples_sigma2_err[iter] <- curr_sigma2_err
  samples_sigma2_spaz[iter] <- curr_sigma2_spaz
  samples_phi[iter] <- curr_phi

  setTxtProgressBar(pb, iter)
}
close(pb)
cat("\nTempo Gibbs Sampler NNGP:", round(Sys.time() - start_time, 2), "secondi\n")

# ==============================================================================
# 5. DIAGNOSTICA FINALE E PLOT
# ==============================================================================
mu_post_mean <- colMeans(samples_mu[burn_in:n_iter, ])

par(mfrow = c(2, 2))

# 1. Scatterplot mu
plot_limits <- range(c(mu_true, mu_post_mean))
plot(mu_true, mu_post_mean, pch = 16, col = "royalblue",
     xlab = "Vero mu", ylab = "Stima a posteriori mu",
     main = paste("Recupero Campo Medio (N=", N, ")", sep=""),
     xlim = plot_limits, ylim = plot_limits)
abline(0, 1, col = "red", lwd = 2, lty = 2)

# 2. Traceplot Varianza Errore
plot(burn_in:n_iter, samples_sigma2_err[burn_in:n_iter], type = "l", col = "grey30",
     xlab = "Iterazione", ylab = "sigma2_err",
     main = paste("Varianza Errore (Vero:", sigma2_errore_true, ")"),
     ylim = c(0.05, 0.2))
abline(h = sigma2_errore_true, col = "red", lwd = 2)

# 3. Traceplot Varianza Spaziale
plot(burn_in:n_iter, samples_sigma2_spaz[burn_in:n_iter], type = "l", col = "darkorange",
     xlab = "Iterazione", ylab = "sigma2_spaz",
     main = paste("Varianza Spaziale (Vero:", sigma2_spaziale_true, ")"),
     ylim = c(1.0, 3.5))
abline(h = sigma2_spaziale_true, col = "red", lwd = 2)

# 4. Traceplot Decadimento Spaziale
plot(burn_in:n_iter, samples_phi[burn_in:n_iter], type = "l", col = "cyan4",
     xlab = "Iterazione", ylab = "phi",
     main = paste("Decadimento Spaziale (Vero:", phi_true, ")"),
     ylim = c(0.5, 2.0))
abline(h = phi_true, col = "red", lwd = 2)

par(mfrow = c(1, 1))

# --- Riassunto Numerico A Posteriori ---
cat("\n=== RIASSUNTO POSTERIORE (Media) ===\n")
cat(sprintf("sigma2_spaziale -> Stimato: %.3f | Vero: %.3f\n", mean(samples_sigma2_spaz[burn_in:n_iter]), sigma2_spaziale_true))
cat(sprintf("phi             -> Stimato: %.3f | Vero: %.3f\n", mean(samples_phi[burn_in:n_iter]), phi_true))
cat(sprintf("sigma2_errore   -> Stimato: %.3f | Vero: %.3f\n", mean(samples_sigma2_err[burn_in:n_iter]), sigma2_errore_true))


# ==============================================================================
# 6. SOLUZIONE ALL'AUTOCORRELAZIONE: THINNING E DIAGNOSTICA ACF
# ==============================================================================

# Definiamo il passo di sfoltimento (teniamo 1 campione ogni 10)
thin_step <- 10
idx_thin <- seq(burn_in, n_iter, by = thin_step)

# Sfoltiamo le catene
samples_sigma2_spaz_thin <- samples_sigma2_spaz[idx_thin]
samples_phi_thin <- samples_phi[idx_thin]
samples_sigma2_err_thin <- samples_sigma2_err[idx_thin]

par(mfrow = c(2, 2))

# --- 1. ACF (Autocorrelation Function) PRIMA del Thinning ---
# Mostra la caduta lenta (la causa delle "onde")
acf(samples_phi[burn_in:n_iter], main = "ACF (phi) - Catena Grezza\n(Elevata Autocorrelazione)")

# --- 2. ACF DOPO il Thinning ---
# Le barre dovrebbero crollare a zero subito, indicando campioni indipendenti
acf(samples_phi_thin, main = paste("ACF (phi) - Dopo Thinning (1 ogni", thin_step, ")"))

# --- 3. Traceplot Sfoltito Varianza Spaziale ---
# Questo è il "bruco denso" che volevamo vedere
plot(idx_thin, samples_sigma2_spaz_thin, type = "l", col = "darkorange",
     xlab = "Iterazione", ylab = "sigma2_spaz",
     main = "Traceplot Sfoltito (sigma2_spaz)",
     ylim = c(1.0, 3.5))
abline(h = sigma2_spaziale_true, col = "red", lwd = 2)

# --- 4. Traceplot Sfoltito Decadimento Spaziale ---
plot(idx_thin, samples_phi_thin, type = "l", col = "cyan4",
     xlab = "Iterazione", ylab = "phi",
     main = "Traceplot Sfoltito (phi)",
     ylim = c(0.5, 2.0))
abline(h = phi_true, col = "red", lwd = 2)

par(mfrow = c(1, 1))

cat("\n=== EFFETTO DEL THINNING ===\n")
cat("Campioni originali post burn-in:", length(burn_in:n_iter), "\n")
cat("Campioni sfoltiti indipendenti utili:", length(idx_thin), "\n")

