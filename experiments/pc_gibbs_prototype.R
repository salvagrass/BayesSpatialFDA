# experiments/pc_gibbs/gibbs_prototype.R
# Main runner for the Marginalised Adaptive Anisotropic SPDE Gibbs Sampler

devtools::load_all()
library(coda)

source("experiments/pc_gibbs/update_mcmc.R")
source("experiments/pc_gibbs/run_sampler.R")

# 1. Setup Mesh & Observations ----
set.seed(42)

dom_length <- 5
domain_bnds <- rbind(
  c(0, 0),
  c(0, dom_length),
  c(dom_length, 0),
  c(dom_length, dom_length)
)

# Parameters studied for approx 500 vertices
mesh <- fm_mesh_2d(
  loc.domain = domain_bnds,
  max.edge = c(0.35, 1.2),
  offset = c(0, 0.08),
  cutoff = 0.25
)

V <- mesh$n
M <- 10L

# Observations not distributed uniformly
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

plot(mesh, main = "Mesh and observations")
points(
  loc_obs,
  pch = 16,
  cex = 0.25,
  col = adjustcolor("#D95F02", alpha.f = 0.25)
)


# 2. Generate sintetic data

matrices <- build_spde_matrices(mesh, manifold = "planar")

truth <- list(
  kappa = 2.0,
  sigma_beta = 0.6,
  v1 = 0.6,
  v2 = 0.4,
  sigma_e = 0.15
)


# Setting up the parameters for the mean field, we generate a plausible function for mu
mu_matrices <- build_spde_matrices(mesh, manifold = "planar")
Q_mu_sim <- build_precision_matrix(
  # Data for Q_mu_sim are chosen with the phenomenon in mind
  spde_matrices = mu_matrices,
  kappa = 0.55,
  sigma_beta = 4.5,
  v1 = 0.7,
  v2 = 0.0
)

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

# 3. Priors & Initial State (deliberately offset from truth) ----
priors <- list(
  lambda_kappa = 0.5,
  lambda_sigma_beta = 0.5,
  sigma_v = 1.0,
  sigma2_epsilon = list(a = 2.0, b = 0.1),
  mu = list(mean = 15.0, sigma2_mu = 200)
)

init_state <- list(
  theta = list(
    kappa = 1.0,
    sigma_beta = 2,
    v1 = 0.0,
    v2 = 0.0
  ),
  sigma2_epsilon = 0.2,
  mu = rep(15.0, V),
  beta_mat = matrix(0.0, nrow = V, ncol = M)
)

# 4. Run Adaptive MCMC Sampler ----
num_iteration <- 5000L
num_burn_in <- 1500L
num_adaptive <- floor(0.3 * num_burn_in)
standard_step <- (0.05^2) * diag(4)

t_start <- proc.time()
fit <- run_gibbs_sampler(
  N = N,
  V = V,
  M = M,
  Phi = sim_data$Phi,
  Y_mat = Y_mat,
  init_state = init_state,
  priors = priors,
  matrices = matrices,
  num_iteration = num_iteration,
  num_burn_in = num_burn_in,
  standard_step = standard_step,
  verbose = TRUE
)
elapsed <- (proc.time() - t_start)[["elapsed"]]

cat(sprintf(
  "Completed %d iterations in %.2f s | Post-burn-in MH acceptance rate: %.1f%%\n",
  num_iteration,
  elapsed,
  100 * fit$acceptance_rate
))

# 5. Trace Plots (Full Trajectory Including Burn-in) ----
# Assumes theta_samples (4 x num_iteration) and sigma2_samples (num_iteration)
# are stored across all iterations; if stored only post-burn-in, set start = num_burn_in + 1.
full_draws <- cbind(
  kappa = fit$theta_samples["kappa", ],
  sigma_beta = fit$theta_samples["sigma_beta", ],
  v1 = fit$theta_samples["v1", ],
  v2 = fit$theta_samples["v2", ],
  sigma2_e = fit$sigma2_samples
)
truth_vec <- c(
  kappa = truth$kappa,
  sigma_beta = truth$sigma_beta,
  v1 = truth$v1,
  v2 = truth$v2,
  sigma2_e = truth$sigma2_e
)

par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
for (param in colnames(full_draws)) {
  plot(
    full_draws[, param],
    type = "l",
    col = "grey30",
    xlab = "Iteration",
    ylab = param,
    main = paste("Trace:", param)
  )
  abline(h = truth_vec[param], col = "firebrick", lwd = 2)
  abline(v = c(num_adaptive, num_burn_in), col = "steelblue", lty = 2)
}
par(mfrow = c(1, 1))

# 6. Post-Burn-In Diagnostics via coda ----
post_idx <- (num_burn_in + 1L):num_iteration
mcmc_post <- mcmc(full_draws[post_idx, ], start = num_burn_in + 1L)

print(summary(mcmc_post))
cat("\nEffective Sample Size (ESS):\n")
print(round(effectiveSize(mcmc_post), 1))

# coda trace + posterior density plots (post-burn-in only)
plot(mcmc_post)

# 7. Spatial Field Recovery Check (mu) ----
mu_post_mean <- rowMeans(fit$mu_samples)
cat(sprintf(
  "\nRMSE on spatial mean field mu: %.4f (prior SD: %.2f, noise SD: %.2f)\n",
  sqrt(mean((mu_post_mean - sim_data$mu_true)^2)),
  sqrt(priors$mu$sigma2_mu),
  truth$sigma_e
))
