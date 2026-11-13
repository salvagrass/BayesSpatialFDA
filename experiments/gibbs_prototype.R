# experiments/01_gibbs_prototype.R
# Main runner for the Marginalised Adaptive Anisotropic SPDE Gibbs Sampler

devtools::load_all(".")
library(coda)

source("experiments/gibbs_sampler/mcmc_updates.R")
source("experiments/gibbs_sampler/run_sampler.R")
source("experiments/gibbs_sampler/diagnostics.R")

# 1. Setup & Data Generation ----
set.seed(42)
x <- seq(0, 5, length.out = 20)
y <- seq(0, 5, length.out = 20)
locs <- as.matrix(expand.grid(x = x, y = y))

mesh <- fm_mesh_2d(loc = locs)
N <- mesh$n
M <- 10L
cat("Nr. of vertices (N):", N, "| Years (M):", M, "\n")

matrices <- build_spde_matrices(mesh, manifold = "planar")

true_tau <- 1.0
true_kappa <- 1.0
true_gamma <- 2.0
true_theta <- pi / 4
true_sigma_e <- sqrt(0.25)

sim <- generate_toy_data(
  mesh = mesh,
  num_year = M,
  mu_prior_mean = 0.5,
  Q_mu = Matrix::Diagonal(N, 1),
  tau = true_tau,
  k = true_kappa,
  gamma = true_gamma,
  theta = true_theta,
  sigma_e = true_sigma_e
)

Y_mat <- sim$y_mat
true_anomaly_mat <- sim$beta_mat - sim$mu_true

truth <- list(
  sigma2 = true_sigma_e^2,
  tau2 = true_tau^2,
  kappa2 = true_kappa^2,
  gamma = true_gamma,
  theta = true_theta,
  beta1 = true_anomaly_mat[1L, 1L]
)

# 2. Run Marginalised Gibbs Sampler ----
res <- run_gibbs_sampler(
  Y_mat = Y_mat,
  matrices = matrices,
  control = list(
    G = 5000L,
    burnin = 1000L,
    thin = 4L,
    adapt_start = 200L,
    adapt_every = 10L,
    init_sd = 0.05,
    scale_shrink = 0.7,
    progress = TRUE
  )
)

# 3. Diagnostics & Validation ----
diag_out <- run_all_diagnostics(
  res = res,
  mesh = mesh,
  Y_mat = Y_mat,
  true_anomaly_mat = true_anomaly_mat,
  truth = truth,
  year = 1L
)
