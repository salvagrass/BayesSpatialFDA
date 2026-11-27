# experiments/pc_gibbs/test_rcpp_sampler.R
# Main runner per il Gibbs Sampler SPDE chiamando il backend C++ via Rcpp

devtools::load_all()
library(coda)
library(fmesher)
library(Matrix)

# 1. SETUP MESH & OSSERVAZIONI ------------------------------------------------
set.seed(42)

dom_length <- 10
domain_bnds <- rbind(
  c(0, 0),
  c(0, dom_length),
  c(dom_length, 0),
  c(dom_length, dom_length)
)

# Mesh su [0, 10]^2 con buffer esterno largo (2.5 unità > 1.5 * rho)
# Mantiene V intorno a ~450-550 vertici totali
mesh <- fm_mesh_2d(
  loc.domain = domain_bnds,
  max.edge = c(0.65, 2.2),
  offset = c(0.3, 2.5),
  cutoff = 0.4
)

V <- mesh$n
M <- 25L

# 10000 osservazioni spaziali non uniformi per ciascuno dei 25 anni
N <- 10000L
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
matrices <- build_spde_matrices(mesh, manifold = "planar")

# Parametri fisici ben condizionati:
# rho = sqrt(8)/1.8 ~ 1.57 (~16% del lato L=10)
# anisotropia a = exp(sqrt(0.5^2 + 0.4^2)) ~ 1.90
truth <- list(
  kappa = 1.8,
  sigma_beta = 1.0,
  v1 = 0.5,
  v2 = 0.4,
  sigma_e = 0.25
)


Q_mu <- fm_matern_precision(mesh, rho = 5, sigma = 3, alpha = 2)

sim_data <- generate_toy_data(
  mesh = mesh,
  num_year = M,
  mu_prior_mean = 15.0,
  Q_mu = 1,
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

Phi <- as(Phi, "dgCMatrix")
Y_mat <- as.matrix(Y_mat)
Gx <- as(matrices$K11, "dgCMatrix")
Gy <- as(matrices$K22, "dgCMatrix")
Gxy <- as(matrices$K12, "dgCMatrix")
C_vec <- as.numeric(Matrix::diag(matrices$C))

cat("\n[R] Starting C++ Gibbs Engine...\n")
t_start <- proc.time()

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
# 4. ANALISI E PLOT CON CODA E GROUND TRUTH (MULTI-CHAIN)
# -----------------------------------------------------------------------------

num_chains <- length(fit)

# Estrazione catene scalari e vettoriali base
sigma_eps_chains <- lapply(fit, function(chain) {
  mcmc(unlist(chain$sigma_epsilon), start = num_burn_in + 1)
})
sigma_eps_mcmc <- mcmc.list(sigma_eps_chains)

theta_chains <- lapply(fit, function(chain) {
  mat <- do.call(rbind, chain$theta)
  colnames(mat) <- c("kappa", "sigma_beta", "v1", "v2")
  mcmc(mat, start = num_burn_in + 1)
})
theta_mcmc <- mcmc.list(theta_chains)

# Estrazione catene per mu[200] e beta[200, 1] (solo per tabella riassuntiva)
mu_200_chains <- lapply(fit, function(chain) {
  mat <- do.call(rbind, chain$mu)
  mcmc(mat[, 200], start = num_burn_in + 1)
})
mu_200_mcmc <- mcmc.list(mu_200_chains)

beta_200_1_chains <- lapply(fit, function(chain) {
  vec <- sapply(chain$beta, function(m) m[200, 1])
  mcmc(vec, start = num_burn_in + 1)
})
beta_200_1_mcmc <- mcmc.list(beta_200_1_chains)

true_sigma2 <- truth$sigma_e^2
true_theta <- c(truth$kappa, truth$sigma_beta, truth$v1, truth$v2)
param_names <- c("kappa", "sigma_beta", "v1", "v2")

# --- 1. Plot Iperparametri ---
cat("\nGenerazione grafici per iperparametri...\n")
par(mfrow = c(4, 2), mar = c(3, 4, 3, 1))

traceplot(sigma_eps_mcmc, main = expression("Traccia - " * sigma[epsilon]^2))
abline(h = true_sigma2, col = "red", lwd = 2, lty = 2)
densplot(sigma_eps_mcmc, main = expression("Densità - " * sigma[epsilon]^2))
abline(v = true_sigma2, col = "red", lwd = 2, lty = 2)

for (i in 1:4) {
  traceplot(
    theta_mcmc[, i, drop = FALSE],
    main = paste("Traccia -", param_names[i])
  )
  abline(h = true_theta[i], col = "red", lwd = 2, lty = 2)
  densplot(
    theta_mcmc[, i, drop = FALSE],
    main = paste("Densità -", param_names[i])
  )
  abline(v = true_theta[i], col = "red", lwd = 2, lty = 2)
}
par(mfrow = c(1, 1))

# --- 2. Tabella Riassuntiva (Mean, Vero, 95% CI) ---
cat("\n======================================================\n")
cat("      CONFRONTO STIME VS GROUND TRUTH\n")
cat("======================================================\n")

# Arrays per calcoli quantilici
se_arr <- as.matrix(sigma_eps_mcmc)
th_arr <- as.matrix(theta_mcmc)
mu_arr <- as.matrix(mu_200_mcmc)
be_arr <- as.matrix(beta_200_1_mcmc)

true_beta_mat <- sim_data$beta_mat

summary_table <- data.frame(
  Parameter = c(
    "sigma2_epsilon",
    "kappa",
    "sigma_beta",
    "v1",
    "v2",
    "mu[200]",
    "beta[200,1]"
  ),
  True_Value = c(
    true_sigma2,
    truth$kappa,
    truth$sigma_beta,
    truth$v1,
    truth$v2,
    sim_data$mu_true[1500],
    true_beta_mat[1500, 1]
  ),
  Mean = unname(c(mean(se_arr), colMeans(th_arr), mean(mu_arr), mean(be_arr))),
  CI_2.5 = unname(c(
    quantile(se_arr, 0.025),
    apply(th_arr, 2, quantile, probs = 0.025),
    quantile(mu_arr, 0.025),
    quantile(be_arr, 0.025)
  )),
  CI_97.5 = unname(c(
    quantile(se_arr, 0.975),
    apply(th_arr, 2, quantile, probs = 0.975),
    quantile(mu_arr, 0.975),
    quantile(be_arr, 0.975)
  ))
)
print(summary_table, row.names = FALSE, digits = 4)

# -----------------------------------------------------------------------------
# 5. CONFRONTO SPAZIALE (Ground Truth vs Posterior)
# -----------------------------------------------------------------------------
cat("\n======================================================\n")
cat("      RICOSTRUZIONE CAMPI SPAZIALI\n")
cat("======================================================\n")

# Estrazione e medie MU
mu_chains_list <- lapply(fit, function(chain) do.call(rbind, chain$mu))
mu_samples <- do.call(rbind, mu_chains_list)
mu_mean <- colMeans(mu_samples)

# Estrazione e medie BETA (Anno 1)
all_beta_matrices <- unlist(lapply(fit, function(x) x$beta), recursive = FALSE)
beta_mean_matrix <- Reduce("+", all_beta_matrices) / length(all_beta_matrices)
beta_year1_mean <- beta_mean_matrix[, 1]

# Centratura campi MU
true_mu_centered <- sim_data$mu_true - mean(sim_data$mu_true)
mu_mean_centered <- mu_mean - mean(mu_mean)

mae_mu <- mean(abs(mu_mean_centered - true_mu_centered))
cat(
  "MU -> MAE:",
  mae_mu,
  "| RMSE:",
  sqrt(mean((mu_mean_centered - true_mu_centered)^2)),
  "\n"
)

mae_beta <- mean(abs(beta_year1_mean - true_beta_mat[, 1]))
cat(
  "BETA (Anno 1) -> MAE:",
  mae_beta,
  "| RMSE:",
  sqrt(mean((beta_year1_mean - true_beta_mat[, 1])^2)),
  "\n"
)

# Plot MU
global_zlim_mu <- range(c(true_mu_centered, mu_mean_centered), na.rm = TRUE)
par(mfrow = c(1, 2))
plot_spde_field(
  mesh,
  true_mu_centered,
  "Vero MU (Centrato)",
  zlim = global_zlim_mu
)
plot_spde_field(
  mesh,
  mu_mean_centered,
  "Stima MU (Centrato)",
  zlim = global_zlim_mu
)

# Plot BETA Anno 1
global_zlim_beta <- range(c(true_beta_mat[, 1], beta_year1_mean), na.rm = TRUE)
par(mfrow = c(1, 2))
plot_spde_field(
  mesh,
  true_beta_mat[, 1],
  "Vero BETA (Anno 1)",
  zlim = global_zlim_beta
)
plot_spde_field(
  mesh,
  beta_year1_mean,
  "Stima BETA (Anno 1)",
  zlim = global_zlim_beta
)

par(mfrow = c(1, 1))
