# ==============================================================================
# PROGETTO BAYESIANO: BASELINE CON MESH REGOLARE, INTERCETTA, E ANOMALIE BETA_i
# ==============================================================================

library(fmesher)
library(MCMCpack)
library(Matrix)
library(devtools)

# --- 1. IMPORTAZIONE DIPENDENZE ---
load_all("R")

# ==============================================================================
# PARTE 1: VALIDAZIONE STATISTICA SU MESH 5x5 (25 VERTICI)
# ==============================================================================
cat("\n--- PARTE 1: ESECUZIONE GIBBS SAMPLER CON BETA_I SU MESH 25 VERTICI ---\n")

set.seed(123)
n_side <- 5
locs_25 <- as.matrix(expand.grid(
  x = seq(0, 5, length.out = n_side),
  y = seq(0, 5, length.out = n_side)
))

toy_mesh_sq <- fm_rcdt_2d_inla(loc = locs_25)
nv <- toy_mesh_sq$n

fem_matrices <- build_spde_matrices(toy_mesh_sq, manifold = "planar")

# Generazione Toy Dataset (25 Anni)
toy_data <- generate_toy_data(
  mesh = toy_mesh_sq,
  sigma_e = 0.1,
  seed = 42
)
Y_mat <- toy_data$y_mat
M <- ncol(Y_mat)
mu_real <- as.numeric(toy_data$mu_true)
a_real <- as.numeric(toy_data$a_true)
beta_real_mat <- toy_data$beta_tilde_mat # Le vere anomalie per i 25 anni

# CENTRATURA DELLA GROUND TRUTH
mean_mu_real <- mean(mu_real)
mu_real <- mu_real - mean_mu_real
a_real <- a_real + mean_mu_real

# Setup MCMC
nu_prior <- nv + 2
Psi_prior <- diag(1, nv)
var_a_prior <- 10.0^2
n_iter <- 5000
burn_in <- 500

samples_mu <- matrix(0, nrow = n_iter, ncol = nv)
samples_a1 <- numeric(n_iter)
samples_beta_nodo1_anno1 <- numeric(n_iter)
samples_Sigma <- vector("list", n_iter)

current_mu <- rep(0, nv)
current_a <- rep(0, M)
current_beta <- matrix(0, nrow = nv, ncol = M)
current_Sigma <- diag(1, nv)
ones <- rep(1, nv)
sigma2_e <- 0.1^2 # Per la baseline assumiamo varianza dell'errore nota per semplicità

# Esecuzione Gibbs Sampler
start_time <- Sys.time()
for (iter in 1:n_iter) {

  Sigma_inv <- solve(current_Sigma)

  # A.1) Intercette (a_i)
  q_a <- (1 / var_a_prior) + as.numeric(t(ones) %*% Sigma_inv %*% ones)
  V_a <- 1 / q_a
  for (i in 1:M) {
    # Il residuo deve ora tenere conto anche dei beta correnti
    resid_i <- Y_mat[, i] - current_mu - current_beta[, i]
    m_a <- V_a * as.numeric(t(ones) %*% Sigma_inv %*% resid_i)
    current_a[i] <- rnorm(1, mean = m_a, sd = sqrt(V_a))
  }

  # A.2) Campo Latente (mu)
  sum_resid <- rep(0, nv)
  for (i in 1:M) {
    sum_resid <- sum_resid + (Y_mat[, i] - current_a[i] - current_beta[, i])
  }
  V_mu_post <- solve(M * Sigma_inv + diag(1e-6, nv))
  m_mu_post <- V_mu_post %*% (Sigma_inv %*% sum_resid)
  current_mu <- MASS::mvrnorm(1, mu = as.vector(m_mu_post), Sigma = V_mu_post)

  # CENTRATURA PER IDENTIFICABILITÀ FRA a_i E mu
  mean_mu <- mean(current_mu)
  current_mu <- current_mu - mean_mu
  current_a <- current_a + mean_mu

  # A.3) AGGIUNTA RIGOROSA: Anomalie Spaziotemporali (beta_i)
  # beta_i ha come prior N(0, current_Sigma) e la likelihood dei dati impone precisione I/sigma2_e
  V_beta_post <- solve(Sigma_inv + diag(1/sigma2_e, nv))
  for(i in 1:M) {
    resid_beta_i <- Y_mat[, i] - current_a[i] - current_mu
    m_beta_post <- V_beta_post %*% (resid_beta_i / sigma2_e)
    current_beta[, i] <- MASS::mvrnorm(1, mu = as.vector(m_beta_post), Sigma = V_beta_post)
  }

  # Centratura dei beta per renderli vere anomalie a media zero
  beta_means <- colMeans(current_beta)
  current_beta <- current_beta - rep(beta_means, each = nv)
  current_a <- current_a + beta_means

  # B) Covarianza (Sigma)
  nu_post <- nu_prior + M
  S_post <- Psi_prior
  for (i in 1:M) {
    # Ora la matrice di covarianza cattura la struttura spaziale reale tramite i beta_i
    S_post <- S_post + (current_beta[, i] %*% t(current_beta[, i]))
  }
  current_Sigma <- MCMCpack::riwish(v = nu_post, S = S_post)

  # Salvataggio
  samples_a1[iter] <- current_a[1]
  samples_mu[iter, ] <- current_mu
  samples_Sigma[[iter]] <- current_Sigma
  samples_beta_nodo1_anno1[iter] <- current_beta[1, 1]

  if (iter %% 1000 == 0) cat("Iterazione:", iter, "/", n_iter, "\n")
}
cat("Tempo Gibbs Sampler:", round(Sys.time() - start_time, 2), "secondi\n")

# --- PLOT DIAGNOSTICI ---
mu_post_mean <- colMeans(samples_mu[burn_in:n_iter, ])
par(mfrow = c(2, 2))

plot_limits <- range(c(mu_real, mu_post_mean))
plot(mu_real, mu_post_mean, pch = 16, col = "royalblue",
     xlab = "Vero mu", ylab = "Stima a posteriori",
     main = "Recupero Campo (Mesh 5x5)", xlim = plot_limits, ylim = plot_limits)
abline(0, 1, col = "red", lwd = 2, lty = 2)

rmse_chain <- sapply(burn_in:n_iter, function(i) sqrt(mean((samples_mu[i, ] - mu_real)^2)))
plot(burn_in:n_iter, rmse_chain, type = "l", col = "darkgreen",
     xlab = "Iterazione", ylab = "RMSE Globale di mu", main = "Convergenza Errore")

var_nodo1 <- sapply(samples_Sigma[burn_in:n_iter], function(mat) mat[1, 1])
plot(burn_in:n_iter, var_nodo1, type = "l", col = "darkorange",
     xlab = "Iterazione", ylab = "Varianza Sigma[1,1]", main = "Traceplot Varianza")

# Sostituito il traceplot dell'intercetta con quello dell'anomalia beta (Nodo 1, Anno 1)
plot(burn_in:n_iter, samples_beta_nodo1_anno1[burn_in:n_iter], type = "l", col = "purple",
     xlab = "Iterazione", ylab = "Beta (Nodo 1, Anno 1)", main = "Traceplot Anomalia Beta")
abline(h = beta_real_mat[1, 1], col = "red", lwd = 2)
par(mfrow = c(1, 1))

# ==============================================================================
# PARTE 2: ESPLOSIONE COMPUTAZIONALE SU GRIGLIE SCALARI (O(N^3))
# ==============================================================================
cat("\n--- PARTE 2: STRESS TEST COMPUTAZIONALE (GRIGLIE REGOLARI) ---\n")

grid_sides <- c(5, 10, 20, 30, 40, 50)
mesh_sizes <- grid_sides^2
inversion_times <- numeric(length(mesh_sizes))
n_rep <- 5

for (i in seq_along(mesh_sizes)) {
  N <- mesh_sizes[i]
  times_rep <- numeric(n_rep)

  for(r in 1:n_rep) {
    dummy_mat <- crossprod(matrix(rnorm(N * N), N, N)) + diag(0.1, N)

    t_start <- proc.time()["elapsed"]
    dummy_inv <- solve(dummy_mat)
    t_end <- proc.time()["elapsed"]

    times_rep[r] <- as.numeric(t_end - t_start)
  }

  inversion_times[i] <- median(times_rep)
  cat("Griglia", grid_sides[i], "x", grid_sides[i], "(Vertici:", N, ") - Tempo mediano:", round(inversion_times[i], 4), "sec\n")
}

# --- PLOT SCALABILITÀ ---
plot(mesh_sizes, inversion_times, type = "b", pch = 16, col = "darkred", lwd = 2,
     main = "Esplosione Computazionale su Mesh Regolari",
     xlab = "Numero di Vertici della Mesh (N)", ylab = "Tempo mediano inversione (sec)")

N_last <- mesh_sizes[length(mesh_sizes)]
T_last <- inversion_times[length(mesh_sizes)]
c_2 <- T_last / (N_last^2)
c_3 <- T_last / (N_last^3)

lines(mesh_sizes, c_2 * mesh_sizes^2, col = "forestgreen", lty = 2, lwd = 2)
lines(mesh_sizes, c_3 * mesh_sizes^3, col = "blue", lty = 2, lwd = 2)

legend("topleft",
       legend = c("Tempo empirico (Mediana)",
                  expression("Crescita Teorica " ~ O(N^2)),
                  expression("Crescita Teorica " ~ O(N^3))),
       col = c("darkred", "forestgreen", "blue"),
       lwd = 2, lty = c(1, 2, 2), pch = c(16, NA, NA), bty = "n")

