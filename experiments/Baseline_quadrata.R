# ==============================================================================
# TEST DEFINITIVO: ACCURATEZZA E TEMPI CON CURVA O(N^3) ANCORATA
# ==============================================================================

library(fmesher)
library(MCMCpack)
library(Matrix)
library(devtools)

# --- 1. IMPORTAZIONE DIPENDENZE ---
load_all("R")

# --- 2. CONFIGURAZIONE DEL BENCHMARK ---
grid_sizes <- c(10, 15, 20, 30)
M_years <- 25
n_iter <- 2000
burn_in <- 500

real_N_list <- numeric(length(grid_sizes))
time_results <- numeric(length(grid_sizes))

all_mu_real <- numeric(0)
all_mu_post <- numeric(0)

cat("========================================================
")
cat(" ESECUZIONE GIBBS SAMPLER SU 4 MESH QUADRATE
")
cat("========================================================

")

# --- 3. CICLO SUI DIVERSI NUMERI DI NODI ---
for (idx in seq_along(grid_sizes)) {

  gs <- grid_sizes[idx]
  cat("--------------------------------------------------------
")
  cat(sprintf(" RUN %d/4: GRIGLIA %dx%d
", idx, gs, gs))

  # A) Creazione Mesh Quadrata
  x_seq <- seq(0, 1, length.out = gs)
  y_seq <- seq(0, 1, length.out = gs)
  locs_square <- as.matrix(expand.grid(x = x_seq, y = y_seq))
  mesh_current <- fm_rcdt_2d_inla(loc = locs_square)

  N <- mesh_current$n
  real_N_list[idx] <- N
  cat(sprintf(" -> Vertici reali della mesh: %d
", N))

  # B) Generazione Dati
  toy_data <- generate_toy_data(
    mesh = mesh_current, num_year = M_years,
    a_prior_mean = 2.0, a_prior_sd = 10.0, sigma_e = 0.1, seed = 42
  )
  Y_mat <- toy_data$y_mat
  mu_real <- as.numeric(toy_data$mu_true)

  # C) Setup Gibbs Sampler
  nu_prior <- N + 2
  Psi_prior <- diag(1, N)
  var_a_prior <- 10.0^2

  samples_mu <- matrix(0, nrow = n_iter, ncol = N)

  current_mu <- rep(0, N)
  current_a <- rep(0, M_years)
  current_Sigma <- diag(1, N)
  ones <- rep(1, N)

  # D) Esecuzione Gibbs
  t_start <- Sys.time()

  for (iter in 1:n_iter) {
    Sigma_inv <- solve(current_Sigma)

    q_a <- (1 / var_a_prior) + as.numeric(t(ones) %*% Sigma_inv %*% ones)
    V_a <- 1 / q_a
    for (i in 1:M_years) {
      resid_i <- Y_mat[, i] - current_mu
      m_a <- V_a * as.numeric(t(ones) %*% Sigma_inv %*% resid_i)
      current_a[i] <- rnorm(1, mean = m_a, sd = sqrt(V_a))
    }

    sum_resid <- rep(0, N)
    for (i in 1:M_years) {
      sum_resid <- sum_resid + (Y_mat[, i] - current_a[i])
    }
    V_mu_post <- solve(M_years * Sigma_inv + diag(1e-6, N))
    m_mu_post <- V_mu_post %*% (Sigma_inv %*% sum_resid)
    current_mu <- MASS::mvrnorm(1, mu = as.vector(m_mu_post), Sigma = V_mu_post)

    nu_post <- nu_prior + M_years
    S_post <- Psi_prior
    for (i in 1:M_years) {
      resid_vec <- Y_mat[, i] - current_a[i] - current_mu
      S_post <- S_post + (resid_vec %*% t(resid_vec))
    }
    current_Sigma <- MCMCpack::riwish(v = nu_post, S = S_post)

    samples_mu[iter, ] <- current_mu
  }

  t_end <- Sys.time()
  exec_time <- as.numeric(difftime(t_end, t_start, units = "secs"))
  time_results[idx] <- exec_time
  cat(sprintf(" -> Tempo di esecuzione: %.2f secondi

", exec_time))

  mu_post_mean <- colMeans(samples_mu[burn_in:n_iter, ])
  all_mu_real <- c(all_mu_real, mu_real)
  all_mu_post <- c(all_mu_post, mu_post_mean)
}


# ==============================================================================
# 4. PLOT (LAYOUT 1x3 CON ANCORAGGIO ALL'ULTIMO PUNTO)
# ==============================================================================

par(mfrow = c(1, 3))

# --- PLOT 1: Accuratezza Globale (Tutte le mesh combinate) ---
plot_limits <- range(c(all_mu_real, all_mu_post))
plot(all_mu_real, all_mu_post, pch = 16, col = "royalblue",
     xlab = "Vero valore di mu", ylab = "Stima a posteriori",
     main = "Accuratezza Gibbs (25-100 Nodi)",
     xlim = plot_limits, ylim = plot_limits)
abline(0, 1, col = "red", lwd = 2, lty = 2)

# --- CALCOLO CURVA TEORICA ANCORATA AL PUNTO PIÙ GRANDE ---
# Invece di usare time_results[1], usiamo l'ultimo punto per non distorcere
# la curva con i costi fissi (overhead) di R sulle matrici piccole.
n_last <- length(real_N_list)
curva_teorica <- time_results[n_last] * (real_N_list / real_N_list[n_last])^3

# --- PLOT 2: Crescita dei Tempi (Lineare) ---
y_max_lin <- max(c(time_results, curva_teorica))
plot(real_N_list, time_results, type = "b", pch = 16, col = "darkred",
     xlab = "Numero di Vertici della Mesh (N)", ylab = "Tempo Totale MCMC (sec)",
     main = "Crescita dei Tempi (Lineare)", ylim = c(0, y_max_lin))

lines(real_N_list, curva_teorica, col = "blue", lty = 2, lwd = 2)
legend("topleft", legend = c("Tempo reale MCMC", "Crescita teorica O(N^3)"),
       col = c("darkred", "blue"), lty = c(1, 2), pch = c(16, NA), cex = 0.8)

# --- PLOT 3: Crescita dei Tempi (Log-Log) ---
y_limits_log <- range(c(log10(time_results), log10(curva_teorica)))
plot(log10(real_N_list), log10(time_results), type = "b", pch = 16, col = "darkgreen",
     xlab = "Log10(Vertici)", ylab = "Log10(Tempo in sec)",
     main = "Crescita dei Tempi (Log-Log)", ylim = y_limits_log)

lines(log10(real_N_list), log10(curva_teorica), col = "blue", lty = 2, lwd = 2)
legend("topleft", legend = c("Dati empirici", "Pendenza Teorica (3)"),
       col = c("darkgreen", "blue"), lty = c(1, 2), pch = c(16, NA), cex = 0.8)

par(mfrow = c(1, 1))

