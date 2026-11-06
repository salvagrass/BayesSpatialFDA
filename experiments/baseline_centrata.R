# ==============================================================================
# PROGETTO BAYESIANO: BASELINE CON INTERCETTA E COLLASSO INVERSE-WISHART
# ==============================================================================

library(fmesher)
library(MCMCpack)
library(Matrix)
library(devtools)

# --- 1. IMPORTAZIONE DIPENDENZE E COMPILAZIONE ---
load_all("R")

# --- 2. CREAZIONE MESH E INTEGRAZIONE C++ ---
set.seed(123)
locs <- matrix(runif(50), ncol = 2)

toy_mesh <- fm_rcdt_2d_inla(loc = locs)
nv <- toy_mesh$n

fem_matrices <- build_spde_matrices(toy_mesh, manifold = "planar")
M_mass <- fem_matrices$C
K_stiffness <- fem_matrices$K

# --- 3. GENERAZIONE TOY DATASET (25 ANNI) ---
toy_data <- generate_toy_data(
  mesh = toy_mesh,
  num_year = 25,
  a_prior_mean = 2.0,
  a_prior_sd = 10.0,
  sigma_e = 0.1,
  seed = 42
)

Y_mat <- toy_data$y_mat # Ora è una matrice 41 x 25
M <- ncol(Y_mat) # M = 25

mu_real <- as.numeric(toy_data$mu_true)
a_real <- as.numeric(toy_data$a_true)

# --- 4. CONFIGURAZIONE PRIORS E SETUP MCMC ---
nu_prior <- nv + 2
Psi_prior <- diag(1, nv)

# Prior per l'intercetta globale (Normale con media 0 e sd 10)
var_a_prior <- 10.0^2

n_iter <- 20000
burn_in <- 5000

samples_mu <- matrix(0, nrow = n_iter, ncol = nv)
samples_a1 <- numeric(n_iter) # Salviamo l'intercetta del primo anno per test
samples_Sigma <- vector("list", n_iter)

# Inizializzazioni
current_mu <- rep(0, nv)
current_a <- rep(0, M) # Vettore di 25 intercette
current_Sigma <- diag(1, nv)
ones <- rep(1, nv)

# --- 5. GIBBS SAMPLER PER M ANNI (LA BASELINE COMPLETA) ---
cat("\nAvvio Gibbs Sampler su matrice", nv, "x", nv, "con", M, "anni...\n")
start_time <- Sys.time()

for (iter in 1:n_iter) {

  # COLLASSO COMPUTAZIONALE O(N^3): L'inversione della densa
  Sigma_inv <- solve(current_Sigma)

  # A.1) Aggiornamento delle 25 Intercette Globali
  # Calcoliamo la precisione per a_i (somma degli elementi di Sigma_inv)
  q_a <- (1 / var_a_prior) + as.numeric(t(ones) %*% Sigma_inv %*% ones)
  V_a <- 1 / q_a

  for (i in 1:M) {
    resid_i <- Y_mat[, i] - current_mu
    m_a <- V_a * as.numeric(t(ones) %*% Sigma_inv %*% resid_i)
    current_a[i] <- rnorm(1, mean = m_a, sd = sqrt(V_a))
  }

  # A.2) Aggiornamento Campo Latente (mu)
  # Sommiamo i residui spaziali di tutti i 25 anni
  sum_resid <- rep(0, nv)
  for (i in 1:M) {
    sum_resid <- sum_resid + (Y_mat[, i] - current_a[i])
  }

  # La precisione di mu scala con M anni. Aggiungiamo 1e-6 per stabilità numerica.
  V_mu_post <- solve(M * Sigma_inv + diag(1e-6, nv))
  m_mu_post <- V_mu_post %*% (Sigma_inv %*% sum_resid)
  current_mu <- MASS::mvrnorm(1, mu = as.vector(m_mu_post), Sigma = V_mu_post)

  # B) Aggiornamento Covarianza (Sigma) tramite Inverse-Wishart
  nu_post <- nu_prior + M
  S_post <- Psi_prior

  # Sommiamo i prodotti esterni per ogni anno
  for (i in 1:M) {
    resid_vec <- Y_mat[, i] - current_a[i] - current_mu
    S_post <- S_post + (resid_vec %*% t(resid_vec))
  }
  current_Sigma <- MCMCpack::riwish(v = nu_post, S = S_post)

  # Salvataggio
  samples_a1[iter] <- current_a[1]
  samples_mu[iter, ] <- current_mu
  samples_Sigma[[iter]] <- current_Sigma

  if (iter %% 500 == 0) cat("Iterazione:", iter, "/", n_iter, "\n")
}
end_time <- Sys.time()
cat("Tempo di esecuzione:", round(end_time - start_time, 2), "secondi\n\n")

# ==============================================================================
# NUOVA DIAGNOSTICA OGGETTIVA E GLOBALE
# ==============================================================================

# 1. Recupero del campo spaziale globale
mu_post_mean <- colMeans(samples_mu[burn_in:n_iter, ])

par(mfrow = c(2, 2))

# PLOT 1: Scatterplot Globale (Verità vs Stima)
plot_limits <- range(c(mu_real, mu_post_mean))
plot(mu_real, mu_post_mean, pch = 16, col = "royalblue",
     xlab = "Vero valore di mu (Ground Truth)", ylab = "Stima a posteriori (media)",
     main = "Recupero del Campo Spaziale", xlim = plot_limits, ylim = plot_limits)
abline(0, 1, col = "red", lwd = 2, lty = 2)

# PLOT 2: RMSE cumulativo (Errore globale nel tempo)
# Calcola l'errore quadratico medio tra l'intero vettore mu stimato e quello reale ad ogni step
rmse_chain <- sapply(burn_in:n_iter, function(i) {
  sqrt(mean((samples_mu[i, ] - mu_real)^2))
})
plot(burn_in:n_iter, rmse_chain, type = "l", col = "darkgreen",
     xlab = "Iterazione", ylab = "RMSE Globale",
     main = "Convergenza dell'Errore Spaziale")

# PLOT 3: Traceplot di un parametro di covarianza (Sigma[1,1])
var_nodo1 <- sapply(samples_Sigma[burn_in:n_iter], function(mat) mat[1, 1])
plot(burn_in:n_iter, var_nodo1, type = "l", col = "darkorange",
     xlab = "Iterazione", ylab = "Varianza (Sigma[1,1])",
     main = "Traceplot Varianza")

# PLOT 4: Traceplot dell'intercetta globale (a_1)
plot(burn_in:n_iter, samples_a1[burn_in:n_iter], type = "l", col = "purple",
     xlab = "Iterazione", ylab = "Intercetta Anno 1",
     main = "Traceplot Intercetta Temporale")
abline(h = a_real[1], col = "red", lwd = 2)

par(mfrow = c(1, 1))

