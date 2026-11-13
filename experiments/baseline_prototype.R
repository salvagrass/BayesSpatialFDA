# ==============================================================================
# TEST DEFINITIVO: ACCURATEZZA E TEMPI GIBBS SAMPLER OTTIMIZZATO (SPDE/C++)
# ==============================================================================

library(fmesher)
library(Matrix)
library(devtools)

# --- 1. IMPORTAZIONE DIPENDENZE ---
load_all("R")

# --- 2. IMPORTAZIONE FUNZIONI GIBBS OTTIMIZZATE (Dal prototipo dei colleghi) ---
update_a <- function(N, M, mu_a_old, tau_a_old, Y_mat, curr_mu, curr_beta, curr_sigma2) {
  tau_data <- N / curr_sigma2
  tau_a_new <- tau_a_old + tau_data
  residual_mat <- Y_mat - curr_mu - curr_beta
  mean_data <- colMeans(residual_mat)
  mean_a_new <- (tau_a_old * mu_a_old + mean_data * tau_data) / tau_a_new
  rnorm(M, mean_a_new, sqrt(1 / tau_a_new))
}

update_mu <- function(N, M, mean_mu_old, prior_mu_prec, Y_mat, curr_a, curr_beta, curr_sigma2) {
  q_data <- Diagonal(N, M / curr_sigma2)
  mean_data <- 1 / M * rowSums(Y_mat - rep(curr_a, each = N) - curr_beta)
  new_q <- q_data + prior_mu_prec
  rhs <- as.vector(q_data %*% mean_data + prior_mu_prec %*% mean_mu_old)
  chol_factor <- chol(new_q)
  new_mean <- as.vector(solve(chol_factor, solve(t(chol_factor), rhs)))
  sample_normal_vector(new_mean, new_q, is_precision = TRUE)
}

update_beta <- function(N, M, beta_prior_prec, Y_mat, curr_a, curr_mu, curr_sigma2) {
  new_q <- beta_prior_prec + Diagonal(N, 1 / curr_sigma2)
  rhs <- (Y_mat - rep(curr_a, each = N) - curr_mu) / curr_sigma2
  chol_factor <- chol(new_q)
  new_mean <- solve(chol_factor, solve(t(chol_factor), rhs))
  Z_mat <- matrix(rnorm(N * M), nrow = N, ncol = M)
  as.matrix(new_mean + solve(chol_factor, Z_mat))
}

metropolis_hastings <- function(curr_tau2, curr_kappa2, curr_beta, Q_prior_beta, C_mat, K_iso, diag_C, prop_sd = 0.08) {
  M <- ncol(curr_beta)
  log_tau2_prop <- rnorm(1, log(curr_tau2), sd = prop_sd)
  log_kappa2_prop <- rnorm(1, log(curr_kappa2), sd = prop_sd)
  tau2_prop <- exp(log_tau2_prop)
  kappa2_prop <- exp(log_kappa2_prop)
  core_prop <- kappa2_prop * C_mat + K_iso
  Q_prop <- forceSymmetric(tau2_prop * core_prop %*% (Diagonal(x = 1 / diag_C) %*% core_prop))

  chol_curr <- chol(Q_prior_beta)
  chol_prop <- chol(Q_prop)

  log_det_curr <- 2 * sum(log(diag(chol_curr)))
  log_det_prop <- 2 * sum(log(diag(chol_prop)))
  quad_curr <- sum((chol_curr %*% curr_beta)^2)
  quad_prop <- sum((chol_prop %*% curr_beta)^2)

  log_lik_curr <- 0.5 * M * log_det_curr - 0.5 * quad_curr
  log_lik_prop <- 0.5 * M * log_det_prop - 0.5 * quad_prop

  log_prior_curr <- dnorm(log(curr_tau2), 0, 5, log = TRUE) + dnorm(log(curr_kappa2), 0, 5, log = TRUE)
  log_prior_prop <- dnorm(log_tau2_prop, 0, 5, log = TRUE) + dnorm(log_kappa2_prop, 0, 5, log = TRUE)

  accepted <- FALSE
  if (log(runif(1)) < (log_lik_prop + log_prior_prop - log_lik_curr - log_prior_curr)) {
    curr_tau2 <- tau2_prop
    curr_kappa2 <- kappa2_prop
    Q_prior_beta <- Q_prop
    accepted <- TRUE
  }
  list(tau2 = curr_tau2, kappa2 = curr_kappa2, Q_beta = Q_prior_beta, accepted = accepted)
}

# --- 3. CONFIGURAZIONE DEL BENCHMARK ---
# Spingiamo fino a 900 nodi (10x10, 15x15, 20x20, 30x30) per mostrare la potenza
grid_sizes <- c(10, 15, 20, 30)
M_years <- 25
n_iter <- 1000 # Riduco a 1000 iterazioni per un benchmark rapido
burn_in <- 250

real_N_list <- numeric(length(grid_sizes))
time_results <- numeric(length(grid_sizes))

all_mu_real <- numeric(0)
all_mu_post <- numeric(0)

cat("========================================================\n")
cat(" ESECUZIONE GIBBS SAMPLER OTTIMIZZATO (SPDE/SPARSO)\n")
cat("========================================================\n\n")

# --- 4. CICLO SUI DIVERSI NUMERI DI NODI ---
for (idx in seq_along(grid_sizes)) {

  gs <- grid_sizes[idx]
  cat("--------------------------------------------------------\n")
  cat(sprintf(" RUN %d/4: GRIGLIA %dx%d \n", idx, gs, gs))

  # A) Creazione Mesh Quadrata
  x_seq <- seq(0, 1, length.out = gs)
  y_seq <- seq(0, 1, length.out = gs)
  locs_square <- as.matrix(expand.grid(x = x_seq, y = y_seq))
  mesh_current <- fm_rcdt_2d_inla(loc = locs_square)

  N <- mesh_current$n
  real_N_list[idx] <- N
  cat(sprintf(" -> Vertici reali della mesh: %d \n", N))

  # B) Costruzione Matrici Sparse SPDE (Il cuore dell'ottimizzazione)
  matrices <- build_spde_matrices(mesh_current, manifold = "planar")
  K_iso <- matrices$K
  C_mat <- matrices$C
  diag_C <- Matrix::diag(C_mat)

  # C) Generazione Dati
  sim <- generate_toy_data(
    mesh = mesh_current, num_year = M_years,
    a_prior_mean = 5, a_prior_sd = 2.0,
    mu_prior_mean = 15.0, Q_mu = Diagonal(N, 1),
    sigma_e = sqrt(0.25)
  )
  Y_mat <- sim$y_mat
  mu_real <- sim$mu_true

  # D) Setup Gibbs Sampler Ottimizzato
  prior_a_prec <- 1 / 100
  prior_a_mean <- 0
  prior_ig_a <- 0.01
  prior_ig_b <- 0.01
  prior_mu_prec <- Diagonal(N, 1 / 100)
  prior_mu_mean <- rep(0, N)

  samples_mu <- matrix(0, nrow = n_iter, ncol = N)

  curr_a <- rep(0, M_years)
  curr_mu <- rep(0, N)
  curr_beta <- matrix(0, nrow = N, ncol = M_years)
  curr_sigma2 <- 1.0
  curr_tau2 <- 0.5
  curr_kappa2 <- 0.5

  core_mat <- curr_kappa2 * C_mat + K_iso
  Q_prior_beta <- forceSymmetric(curr_tau2 * core_mat %*% (Diagonal(x = 1 / diag_C) %*% core_mat))

  # E) Esecuzione Gibbs Ottimizzato
  t_start <- Sys.time()

  for (iter in 1:n_iter) {
    # Step 1: Intercette
    curr_a <- update_a(N, M_years, prior_a_mean, prior_a_prec, Y_mat, curr_mu, curr_beta, curr_sigma2)
    a_bar <- mean(curr_a)
    curr_a <- curr_a - a_bar
    curr_mu <- curr_mu + a_bar

    # Step 2: Campo Globale (mu) - Usa matrici sparse di precisione!
    curr_mu <- update_mu(N, M_years, prior_mu_mean, prior_mu_prec, Y_mat, curr_a, curr_beta, curr_sigma2)

    # Step 3: Campi Latenti Spazio-Temporali (beta_tilde_i)
    curr_beta <- update_beta(N, M_years, Q_prior_beta, Y_mat, curr_a, curr_mu, curr_sigma2)
    beta_means <- colMeans(curr_beta)
    curr_beta <- curr_beta - rep(beta_means, each = N)
    curr_a <- curr_a + beta_means

    # Step 4: Varianza Errore
    eps_mat <- Y_mat - rep(curr_a, each = N) - curr_mu - curr_beta
    curr_sigma2 <- 1 / rgamma(1, shape = (N * M_years) / 2 + prior_ig_a, rate = sum(eps_mat^2) / 2 + prior_ig_b)

    # Step 5: Iperparametri SPDE (Sostituisce Inverse-Wishart)
    mh_step <- metropolis_hastings(curr_tau2, curr_kappa2, curr_beta, Q_prior_beta, C_mat, K_iso, diag_C, prop_sd = 0.08)
    curr_tau2 <- mh_step$tau2
    curr_kappa2 <- mh_step$kappa2
    Q_prior_beta <- mh_step$Q_beta

    samples_mu[iter, ] <- curr_mu
  }

  t_end <- Sys.time()
  exec_time <- as.numeric(difftime(t_end, t_start, units = "secs"))
  time_results[idx] <- exec_time
  cat(sprintf(" -> Tempo esecuzione OTTIMIZZATO: %.2f secondi\n\n", exec_time))

  # --- Sostituisci la parte finale del ciclo for con questo: ---
  mu_post_mean <- colMeans(samples_mu[burn_in:n_iter, ])

  # CENTRATURA: Rimuoviamo la media globale per confrontare solo la struttura spaziale
  mu_real_centered <- mu_real - mean(mu_real)
  mu_post_centered <- mu_post_mean - mean(mu_post_mean)

  all_mu_real <- c(all_mu_real, mu_real_centered)
  all_mu_post <- c(all_mu_post, mu_post_centered)
}

# ==============================================================================
# 5. PLOT DIAGNOSTICI: L'IMPATTO DELL'APPROCCIO SPARSO (SPDE)
# ==============================================================================

par(mfrow = c(1, 3))

# --- PLOT 1: Accuratezza Globale (SPDE) ---
plot_limits <- range(c(all_mu_real, all_mu_post))
plot(all_mu_real, all_mu_post, pch = 16, col = "royalblue",
     xlab = "Vero valore di mu", ylab = "Stima a posteriori",
     main = "Accuratezza SPDE (100-900 Nodi)",
     xlim = plot_limits, ylim = plot_limits)
abline(0, 1, col = "red", lwd = 2, lty = 2)

# --- CALCOLO CURVE TEORICHE: SPARSO O(N^1.5) vs DENSO O(N^3) ---
n_last <- length(real_N_list)
# Il nuovo modello NON è più O(N^3). Scalerà circa come O(N^1.5)
curva_teorica_sparsa <- time_results[n_last] * (real_N_list / real_N_list[n_last])^1.5
curva_teorica_densa_per_confronto <- time_results[n_last] * (real_N_list / real_N_list[n_last])^3

# --- PLOT 2: Crescita dei Tempi (Lineare) ---
y_max_lin <- max(c(time_results, curva_teorica_densa_per_confronto))
plot(real_N_list, time_results, type = "b", pch = 16, col = "darkgreen",
     xlab = "Numero di Vertici della Mesh (N)", ylab = "Tempo Totale MCMC (sec)",
     main = "Collasso Tempi (SPDE vs Denso)", ylim = c(0, y_max_lin))

lines(real_N_list, curva_teorica_densa_per_confronto, col = "red", lty = 2, lwd = 2) # Proiezione Vecchia Densa
lines(real_N_list, curva_teorica_sparsa, col = "blue", lty = 2, lwd = 2)            # Proiezione Nuova Sparsa
legend("topleft", legend = c("Tempo SPDE reale", "O(N^3) vecchio modello", "O(N^1.5) attesa SPDE"),
       col = c("darkgreen", "red", "blue"), lty = c(1, 2, 2), pch = c(16, NA, NA), cex = 0.8)

# --- PLOT 3: Crescita dei Tempi (Log-Log) ---
y_limits_log <- range(c(log10(time_results), log10(curva_teorica_densa_per_confronto)))
plot(log10(real_N_list), log10(time_results), type = "b", pch = 16, col = "darkgreen",
     xlab = "Log10(Vertici)", ylab = "Log10(Tempo in sec)",
     main = "Scala Log-Log (Appiattimento SPDE)", ylim = y_limits_log)

lines(log10(real_N_list), log10(curva_teorica_densa_per_confronto), col = "red", lty = 2, lwd = 2)
lines(log10(real_N_list), log10(curva_teorica_sparsa), col = "blue", lty = 2, lwd = 2)
legend("topleft", legend = c("Dati empirici SPDE", "Pendenza Densa (3)", "Pendenza Sparsa (1.5)"),
       col = c("darkgreen", "red", "blue"), lty = c(1, 2, 2), pch = c(16, NA, NA), cex = 0.8)

par(mfrow = c(1, 1))

