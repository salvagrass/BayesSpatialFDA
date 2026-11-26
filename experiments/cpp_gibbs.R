# experiments/pc_gibbs/test_rcpp_sampler.R
# Main runner per il Gibbs Sampler SPDE chiamando il backend C++ via Rcpp

devtools::load_all()
library(coda)
library(fmesher)
library(Matrix)

# 1. SETUP MESH & OSSERVAZIONI ------------------------------------------------
set.seed(42)

dom_length <- 5
domain_bnds <- rbind(
  c(0, 0),
  c(0, dom_length),
  c(dom_length, 0),
  c(dom_length, dom_length)
)

# Mesh con ~500 vertici
mesh <- fm_mesh_2d(
  loc.domain = domain_bnds,
  max.edge = c(0.35, 1.2),
  offset = c(0, 0.08),
  cutoff = 0.25
)

V <- mesh$n
M <- 10L

# 10.000 Osservazioni non uniformi
N <- 2000
loc_obs <- cbind(
  x = runif(N, 0, dom_length),
  y = dom_length * rbeta(N, shape1 = 0.65, shape2 = 0.65)
)

cat(
  "Nr. of vertices (V):",
  V,
  "| Years (M):",
  M,
  "| Observations (N): ",
  N,
  "\n"
)

# 2. GENERAZIONE DATI SINTETICI -----------------------------------------------
# Costruiamo le matrici FEM di base dal tuo assembler C++!
matrices <- build_spde_matrices(mesh, manifold = "planar")

truth <- list(
  kappa = 2.0,
  sigma_beta = 0.6,
  v1 = 0.6,
  v2 = 0.4,
  sigma_e = 0.15
)

# Simulazione di mu (Assumiamo tu abbia una funzione build_precision_matrix in R o usi INLA)
mu_matrices <- build_spde_matrices(mesh, manifold = "planar")

Q_mu_sim <- build_precision_matrix(
  spde_matrices = mu_matrices,
  kappa = 0.55,
  sigma_beta = 4.5,
  v1 = 0.7,
  v2 = 0.0
)

# Se per caso la matrice continua a dare noie numeriche per via di arrotondamenti macchina (1e-16),
# puoi forzarne la simmetria in R così:
Q_mu_sim <- Matrix::forceSymmetric(Q_mu_sim)

sim_data <- generate_toy_data(
  mesh = mesh,
  num_year = M,
  mu_prior_mean = 15.0,
  Q_mu = Q_mu_sim,
  kappa = truth$kappa,
  sigma_beta = truth$sigma_beta,
  v1 = truth$v1,
  v2 = truth$v2,
  sigma_e = truth$sigma_e,
  locs = loc_obs,
  seed = 123
)

Y_mat <- sim_data$y_mat
Phi <- sim_data$Phi

# 3. CHIAMATA AL GIBBS SAMPLER C++ (Zero Allocations Architecture) -----------
num_iteration <- 5000L
num_burn_in <- 1500L

# Pulizia dei tipi per RcppEigen
Phi <- as(Phi, "dgCMatrix")
Y_mat <- as.matrix(Y_mat)
Gx <- as(matrices$K11, "dgCMatrix")
Gy <- as(matrices$K22, "dgCMatrix")
Gxy <- as(matrices$K12, "dgCMatrix")
C_vec <- as.numeric(Matrix::diag(matrices$C))

cat("\n[R] Starting C++ Gibbs Engine...\n")
t_start <- proc.time()


# Ora puoi chiamare il tuo campionatore C++
# res <- run_gibbs_sampler(Phi, Y, Gx, Gy, Gxy, ...)

# Passiamo ESATTAMENTE le variabili pulite che abbiamo appena creato!
fit <- .run_gibbs_sampler_cpp(
  Phi = Phi,
  Y = Y_mat,
  Gx = Gx,
  Gy = Gy,
  Gxy = Gxy,
  C = C_vec,
  num_iters = num_iteration,
  num_burn_in = num_burn_in
)

elapsed <- (proc.time() - t_start)[["elapsed"]]
cat(sprintf("\n[R] C++ Sampler Finished in %.2f s\n", elapsed))

# -----------------------------------------------------------------------------
# 4. ANALISI E PLOT CON CODA E GROUND TRUTH
# -----------------------------------------------------------------------------

# --- Preparazione oggetti MCMC ---
sigma_eps_vec <- unlist(fit$sigma_epsilon)
sigma_eps_mcmc <- mcmc(sigma_eps_vec, start = num_burn_in + 1)

theta_mat <- do.call(rbind, fit$theta)
colnames(theta_mat) <- c("kappa", "sigma_beta", "v1", "v2")
theta_mcmc <- mcmc(theta_mat, start = num_burn_in + 1)

# Verità di base
true_sigma2 <- truth$sigma_e^2
true_theta <- c(truth$kappa, truth$sigma_beta, truth$v1, truth$v2)
param_names <- colnames(theta_mat)


# --- 1. Plot Sigma2 Epsilon ---
cat("\nGenerazione grafici per sigma2_epsilon...\n")

# Imposta la finestra grafica su 1 riga e 2 colonne
par(mfrow = c(1, 2))

# Traccia
traceplot(sigma_eps_mcmc, main = expression("Traccia - " * sigma[epsilon]^2))
abline(h = true_sigma2, col = "red", lwd = 2, lty = 2)

# Densità
densplot(sigma_eps_mcmc, main = expression("Densità - " * sigma[epsilon]^2))
abline(v = true_sigma2, col = "red", lwd = 2, lty = 2)


# --- 2. Plot Theta ---
cat("Generazione grafici per il blocco Theta...\n")

# Imposta la finestra grafica su 4 righe e 2 colonne per vedere tutto insieme
# Modifichiamo leggermente i margini (mar) per non schiacciare troppo i titoli
par(mfrow = c(4, 2), mar = c(3, 4, 3, 1))

for (i in 1:4) {
  param <- param_names[i]
  true_val <- true_theta[i]

  # Plot Traccia
  traceplot(theta_mcmc[, i], main = paste("Traccia -", param))
  abline(h = true_val, col = "red", lwd = 2, lty = 2)

  # Plot Densità
  densplot(theta_mcmc[, i], main = paste("Densità -", param))
  abline(v = true_val, col = "red", lwd = 2, lty = 2)
}

# Ripristina la finestra grafica alle impostazioni di default
par(mfrow = c(1, 1))

# --- Sommario Statistiche ---
cat("\n--- Sommario a posteriori (Sigma Epsilon) ---\n")
print(summary(sigma_eps_mcmc))

cat("\n--- Sommario a posteriori (Theta) ---\n")
print(summary(theta_mcmc))
