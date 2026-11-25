# ==============================================================================
# REAL-DATA PLANAR TEST - ANISOTROPIC SPDE PC-GIBBS SAMPLER (MERGED VERSION)
# Model: Y_i = Phi %*% beta_i + eps_i,   beta_i | mu, theta ~ N_V(mu, Q(theta)^-1)
# ==============================================================================

devtools::load_all()
library(Matrix)
library(coda)
library(ggplot2)
library(grid)

source("experiments/pc_gibbs/update_mcmc.R")
source("experiments/pc_gibbs/run_sampler.R")

# ------------------------------------------------------------------------------
# 1. CARICAMENTO DATI REALI E MESH CON BORDO ESTERNO AMPIO (UNITÀ = side_km)
# ------------------------------------------------------------------------------
set.seed(42)

# Mesh con doppio cuscinetto esterno (offset) per allontanare gli effetti di bordo
# di Neumann dal quadrato di interesse; triangoli più larghi fuori per contenere V.
real_data <- prepare_real_planar_data(
  path          = "data/climate_data.RData",
  center        = c(41.9028, 12.4964), # Roma
  shape         = "square",
  side_km       = 2000,
  scale_to_unit = TRUE,
  build_mesh    = TRUE,
  max_edge      = c(0.10, 0.25),       # fitto dentro, più largo nel buffer
  offset        = c(0.05, 0.25),       # bordo esterno ampio (500 km extra per lato)
  cutoff        = 0.04
)

# caso KNN
# real_data <- prepare_real_planar_data(
#   path          = "data/climate_data.RData",
#   center        = c(41.9028, 12.4964), # Roma
#   shape         = "knn",
#   k             = 400L,
#   scale_to_unit = TRUE,
#   build_mesh    = TRUE,
#   max_edge      = c(0.10, 0.25),
#   offset        = c(0.05, 0.25),
#   cutoff        = 0.04
# )

Y_mat    <- real_data$y_mat
Phi      <- real_data$Phi
mesh     <- real_data$mesh
years    <- real_data$years
scale_km <- real_data$meta$scale_factor_km

N <- nrow(Y_mat)
M <- ncol(Y_mat)
V <- mesh$n

empty_idx   <- which(Matrix::colSums(Phi) == 0)
empty_nodes <- length(empty_idx)

cat(sprintf(
  "Dati reali caricati -> N oss: %d | V nodi mesh: %d (%.1f oss/nodo) | M anni: %d (%d-%d) | Nodi vuoti (buffer): %d\n",
  N, V, N / V, M, min(years), max(years), empty_nodes
))
stopifnot(max(abs(Matrix::rowSums(Phi) - 1)) < 1e-10)

coords_km <- real_data$coords * scale_km
bnd_km    <- real_data$boundary * scale_km

# Plot diagnostico della mesh (con bordo rosso e nodi vuoti del buffer)
par(mfrow = c(1, 1))
plot(
  mesh,
  main = sprintf("Planar Mesh & Observations (N = %d, V = %d)", N, V),
  asp = 1
)
polygon(real_data$boundary, border = "red", lwd = 2)
points(real_data$coords, pch = 16, cex = 0.5, col = adjustcolor("steelblue", 0.8))
if (empty_nodes > 0L) {
  points(mesh$loc[empty_idx, 1:2, drop = FALSE], pch = 4, col = "darkorange", lwd = 1.8, cex = 0.85)
}
legend(
  "topright",
  legend = c(
    sprintf("Observations (N = %d)", N),
    "Square boundary",
    if (empty_nodes > 0L) sprintf("Empty nodes (%d)", empty_nodes)
  ),
  col = c("steelblue", "red", if (empty_nodes > 0L) "darkorange"),
  pch = c(16, NA, if (empty_nodes > 0L) 4),
  lty = c(NA, 1, if (empty_nodes > 0L) NA),
  lwd = c(NA, 2, if (empty_nodes > 0L) 2),
  bg = "white", cex = 0.82
)

# ------------------------------------------------------------------------------
# 2. ASSEMBLAGGIO MATRICI SPDE ANISOTROPE (C, K11, K12, K22)
# ------------------------------------------------------------------------------
spde_matrices <- build_spde_matrices(mesh, manifold = "planar")
if (!all(c("C", "K11", "K12", "K22") %in% names(spde_matrices))) {
  stop(
    "build_spde_matrices() ha restituito: ", paste(names(spde_matrices), collapse = ", "),
    ". Attese C, K11, K12, K22 -> controlla R/matrix_assembler.R e rilancia Rcpp::compileAttributes()."
  )
}

# ------------------------------------------------------------------------------
# 3. PC PRIORS FISICHE E INIZIALIZZAZIONE INTELLIGENTE (UNITÀ = scale_km)
# ------------------------------------------------------------------------------
# PC prior su kappa: P(rho < rho0_km) = P(kappa > sqrt(8)/rho0) = alpha_rho
rho0_km <- 200
rho0    <- rho0_km / scale_km
alpha_rho <- 0.05
lambda_kappa <- -(rho0 / sqrt(8)) * log(alpha_rho)

# PC prior su sigma_beta: P(sigma_beta > sig0) = alpha_sig
sig0 <- 5.0
alpha_sig <- 0.05
lambda_sigma_beta <- -log(alpha_sig) / sig0

# Prior su vettore anisotropia v = (v1, v2) ~ N(0, sigma_v^2 I_2)
sigma_v <- 1.0

priors <- list(
  lambda_kappa      = lambda_kappa,
  lambda_sigma_beta = lambda_sigma_beta,
  sigma_v           = sigma_v,
  sigma2_epsilon    = list(a = 2.0, b = 0.1),
  mu                = list(mean = mean(Y_mat), sigma2_mu = 20^2)
)

# Inizializzazione informata via proiezione regolarizzata di Y_bar sulla mesh
# (evita lo shock iniziale su sigma2_epsilon all'iterazione 1)
PhiTPhi_reg <- Matrix::crossprod(Phi) + 1e-3 * Matrix::Diagonal(V)
mu_init     <- as.vector(Matrix::solve(PhiTPhi_reg, Matrix::crossprod(Phi, rowMeans(Y_mat))))
beta_init   <- matrix(mu_init, nrow = V, ncol = M)

init_state <- list(
  theta = list(
    kappa      = sqrt(8) / (1000 / scale_km), # range iniziale = 1000 km
    sigma_beta = 1.0,                         # SD marginale iniziale = 1 °C
    v1         = 0.0,                         # partenza isotropa
    v2         = 0.0
  ),
  sigma2_epsilon = 0.25,
  mu             = mu_init,
  beta_mat       = beta_init
)

# Covarianza iniziale della proposta MH su eta = (log(kappa), log(sigma_beta), v1, v2)
standard_step <- diag(c(0.08, 0.08, 0.10, 0.10)^2)

# ------------------------------------------------------------------------------
# 4. ESECUZIONE ADAPTIVE PC-GIBBS SAMPLER
# ------------------------------------------------------------------------------
num_iteration <- 5000L
num_burn_in   <- 1500L
num_adaptive  <- floor(0.3 * num_burn_in)

cat("\nAvvio Adaptive Anisotropic PC-Gibbs Sampler sui dati reali...\n")
t_start <- proc.time()

# Nota: l'argomento in run_sampler.R su feature/sampler-cpp è al singolare (num_iteration)
fit <- run_gibbs_sampler(
  N             = N,
  V             = V,
  M             = M,
  Phi           = Phi,
  Y_mat         = Y_mat,
  init_state    = init_state,
  priors        = priors,
  matrices      = spde_matrices,
  num_iteration = num_iteration,
  num_burn_in   = num_burn_in,
  standard_step = standard_step,
  verbose       = TRUE
)

elapsed <- (proc.time() - t_start)[["elapsed"]]
cat(sprintf(
  "\nCompletate %d iterazioni in %.2f s | MH acceptance rate (post burn-in): %.1f%%\n",
  num_iteration, elapsed, 100 * fit$acceptance_rate
))

# ------------------------------------------------------------------------------
# 5. DIAGNOSTICA MCMC (CODA, ESS) E RIEPILOGO FISICO A POSTERIORI
# ------------------------------------------------------------------------------
full_draws <- cbind(
  kappa      = fit$theta_samples["kappa", ],
  sigma_beta = fit$theta_samples["sigma_beta", ],
  v1         = fit$theta_samples["v1", ],
  v2         = fit$theta_samples["v2", ],
  sigma2_e   = fit$sigma2_samples
)

n_stored <- nrow(full_draws)
if (n_stored == num_iteration) {
  post_idx   <- (num_burn_in + 1L):num_iteration
  mcmc_start <- num_burn_in + 1L
} else {
  post_idx   <- seq_len(n_stored)
  mcmc_start <- num_burn_in + 1L
}

# Traceplot in scala fisica interpretabile (con ylim tarato sul post-burn-in per sigma_eps)
sig_e_full <- sqrt(full_draws[, "sigma2_e"])
sig_e_ylim <- range(sig_e_full[post_idx]) * c(0.995, 1.3)

par(mfrow = c(2, 3), mar = c(4, 4, 2.2, 1))
phys_traces <- list(
  list(y = sig_e_full,                                        main = "sigma_eps (°C) [zoom]", ylab = "°C",   ylim = sig_e_ylim),
  list(y = (sqrt(8) / full_draws[, "kappa"]) * scale_km,      main = "Range rho (km)",          ylab = "km",   ylim = NULL),
  list(y = full_draws[, "sigma_beta"],                        main = "sigma_beta (°C)",         ylab = "°C",   ylim = NULL),
  list(y = full_draws[, "v1"],                                main = "Anisotropia v1",          ylab = "v1",   ylim = NULL),
  list(y = full_draws[, "v2"],                                main = "Anisotropia v2",          ylab = "v2",   ylim = NULL),
  list(y = sqrt(full_draws[, "v1"]^2 + full_draws[, "v2"]^2), main = "Norma ||v||",             ylab = "||v||", ylim = NULL)
)
for (tr in phys_traces) {
  plot(tr$y, type = "l", col = "grey30", xlab = "Iterazione",
       ylab = tr$ylab, main = tr$main, ylim = tr$ylim)
  if (n_stored == num_iteration) {
    abline(v = c(num_adaptive, num_burn_in), col = c("steelblue", "firebrick"), lty = 2)
  }
}
par(mfrow = c(1, 1))

# Riepilogo coda + ESS
mcmc_post <- coda::mcmc(full_draws[post_idx, ], start = mcmc_start)
print(summary(mcmc_post))
cat("\nEffective Sample Size (ESS):\n")
print(round(coda::effectiveSize(mcmc_post), 1))

# Grandezze fisiche derivate (Slide 10)
theta_post      <- full_draws[post_idx, , drop = FALSE]
kappa_post      <- theta_post[, "kappa"]
sig_b_post      <- theta_post[, "sigma_beta"]
v1_post         <- theta_post[, "v1"]
v2_post         <- theta_post[, "v2"]
sig_e_post      <- sqrt(theta_post[, "sigma2_e"])
rho_km_post     <- (sqrt(8) / kappa_post) * scale_km
v_norm_post     <- sqrt(v1_post^2 + v2_post^2)
aniso_ratio     <- exp(v_norm_post)
aniso_angle_vec <- (0.5 * atan2(mean(v2_post), mean(v1_post))) * (180 / pi)

cat("\n=== RIEPILOGO FISICO A POSTERIORI (MEDIE & IC 95%) ===\n")
cat(sprintf("sigma_eps   (°C)   : mean = %.4f | IC 95%% = [%.4f, %.4f]\n",
            mean(sig_e_post), quantile(sig_e_post, 0.025), quantile(sig_e_post, 0.975)))
cat(sprintf("sigma_beta  (°C)   : mean = %.4f | IC 95%% = [%.4f, %.4f]\n",
            mean(sig_b_post), quantile(sig_b_post, 0.025), quantile(sig_b_post, 0.975)))
cat(sprintf("kappa       (unit) : mean = %.4f | IC 95%% = [%.4f, %.4f]\n",
            mean(kappa_post), quantile(kappa_post, 0.025), quantile(kappa_post, 0.975)))
cat(sprintf("rho         (km)   : mean = %.1f   | IC 95%% = [%.1f, %.1f]\n",
            mean(rho_km_post), quantile(rho_km_post, 0.025), quantile(rho_km_post, 0.975)))
cat(sprintf("||v||              : mean = %.4f | IC 95%% = [%.4f, %.4f]\n",
            mean(v_norm_post), quantile(v_norm_post, 0.025), quantile(v_norm_post, 0.975)))
cat(sprintf("aniso ratio e^||v||: mean = %.3f  | IC 95%% = [%.3f, %.3f]\n",
            mean(aniso_ratio), quantile(aniso_ratio, 0.025), quantile(aniso_ratio, 0.975)))
cat(sprintf("direzione princ.   : %.1f° (0° = Est-Ovest, ±90° = Nord-Sud)\n", aniso_angle_vec))

# ------------------------------------------------------------------------------
# 6. CAMPO MEDIO m(s) E ANOMALIE (GESTIONE DIFENSIVA POST-BURN-IN)
# ------------------------------------------------------------------------------
if (ncol(fit$mu_samples) == num_iteration) {
  mu_post_draws <- fit$mu_samples[, post_idx, drop = FALSE]
} else {
  mu_post_draws <- fit$mu_samples
}
mu_post_mean <- rowMeans(mu_post_draws)
mu_hat_obs   <- as.numeric(Phi %*% mu_post_mean)
y_emp_mean   <- rowMeans(Y_mat)
anom_emp_obs <- sweep(Y_mat, 1, y_emp_mean, "-")

cat(sprintf("\nRMSE campo medio (osservato vs Phi*mu_hat): %.3f °C\n",
            sqrt(mean((y_emp_mean - mu_hat_obs)^2))))

# Identificazione automatica degli anni estremi (più freddo e più caldo nella regione)
mean_emp_by_year <- colMeans(anom_emp_obs)
idx_cold <- which.min(mean_emp_by_year)
idx_warm <- which.max(mean_emp_by_year)
target_idx <- c(cold = idx_cold, warm = idx_warm)

cat(sprintf(
  "Anni estremi selezionati -> Freddo: %d (%+.2f °C) | Caldo: %d (%+.2f °C)\n",
  years[idx_cold], mean_emp_by_year[idx_cold],
  years[idx_warm], mean_emp_by_year[idx_warm]
))

# ------------------------------------------------------------------------------
# 7. PREDIZIONE CONTINUA (80x80) E INCERTEZZA A POSTERIORI (VARIANZA TOTALE)
# ------------------------------------------------------------------------------
n_grid <- 80L
x_seq_unit <- seq(min(real_data$boundary[, 1]), max(real_data$boundary[, 1]), length.out = n_grid)
y_seq_unit <- seq(min(real_data$boundary[, 2]), max(real_data$boundary[, 2]), length.out = n_grid)
grid_unit  <- as.matrix(expand.grid(x = x_seq_unit, y = y_seq_unit))
base_grid  <- expand.grid(x = x_seq_unit * scale_km, y = y_seq_unit * scale_km)

# Valutazione FEM continua (mascheratura automatica compatibile con square e knn)
Phi_grid_all <- as.matrix(fmesher::fm_basis(mesh, loc = grid_unit))

# Taglia sia il buffer esterno sia gli angolini privi di osservazioni
mesh_obs <- fmesher::fm_rcdt_2d_inla(loc = real_data$coords)
in_mesh  <- Matrix::rowSums(fmesher::fm_basis(mesh_obs, loc = grid_unit)) > 0.999

Phi_grid_in  <- Phi_grid_all[in_mesh, , drop = FALSE]
n_in         <- nrow(Phi_grid_in)

n_post  <- nrow(theta_post)
sub_idx <- unique(round(seq(1, n_post, length.out = min(150L, n_post))))
n_sub   <- length(sub_idx)

PtP    <- Matrix::crossprod(Phi)
PtY    <- as.matrix(Matrix::crossprod(Phi, Y_mat))

cond_var_acc  <- numeric(n_in)
mu_grid_draws <- matrix(0, nrow = n_in, ncol = n_sub)
a_cold_draws  <- matrix(0, nrow = n_in, ncol = n_sub)
a_warm_draws  <- matrix(0, nrow = n_in, ncol = n_sub)
f_warm_draws  <- matrix(0, nrow = n_in, ncol = n_sub)
beta_all_acc  <- matrix(0, nrow = V, ncol = M)

for (j in seq_len(n_sub)) {
  s <- sub_idx[j]
  s_mu <- if (ncol(mu_post_draws) == n_post) s else {
    min(ncol(mu_post_draws), max(1L, round(s * ncol(mu_post_draws) / n_post)))
  }
  mu_s <- mu_post_draws[, s_mu]

  Q_th <- build_precision_matrix(
    spde_matrices = spde_matrices,
    kappa         = theta_post[s, "kappa"],
    sigma_beta    = theta_post[s, "sigma_beta"],
    v1            = theta_post[s, "v1"],
    v2            = theta_post[s, "v2"]
  )
  inv_s2 <- 1 / theta_post[s, "sigma2_e"]

  Q_post     <- Matrix::forceSymmetric(inv_s2 * PtP + Q_th)
  Sigma_post <- as.matrix(Matrix::solve(Q_post))

  # Varianza condizionale di phi(s)^T beta_i sui punti interni alla mesh
  cond_var_acc <- cond_var_acc + rowSums((Phi_grid_in %*% Sigma_post) * Phi_grid_in)

  # Medie condizionali di tutte le anomalie (V x M)
  rhs_all <- inv_s2 * sweep(PtY, 1, as.numeric(PtP %*% mu_s), "-")
  m_anom  <- Sigma_post %*% rhs_all
  beta_all_acc <- beta_all_acc + sweep(m_anom, 1, mu_s, "+")

  mu_g <- as.numeric(Phi_grid_in %*% mu_s)
  mu_grid_draws[, j] <- mu_g
  a_cold_draws[, j]  <- as.numeric(Phi_grid_in %*% m_anom[, idx_cold])
  a_warm_draws[, j]  <- as.numeric(Phi_grid_in %*% m_anom[, idx_warm])
  f_warm_draws[, j]  <- mu_g + a_warm_draws[, j]
}

# In-sample RMSE globale (tutti gli anni) sul campo totale
beta_post_mean <- beta_all_acc / n_sub
f_all_obs      <- as.matrix(Phi %*% beta_post_mean)
rmse_all       <- sqrt(mean((Y_mat - f_all_obs)^2))
rmse_2007      <- sqrt(mean((Y_mat[, which(years == 2007L)] - f_all_obs[, which(years == 2007L)])^2))
cat(sprintf("In-sample RMSE campo totale -> tutti gli anni: %.3f °C | anno 2007: %.3f °C\n",
            rmse_all, rmse_2007))

# Ricomposizione vettori completi (NA fuori dalla mesh, utile per KNN)
fill_grid <- function(v_in) {
  out <- rep(NA_real_, nrow(grid_unit))
  out[in_mesh] <- v_in
  out
}
row_var    <- function(mat) rowSums((mat - rowMeans(mat))^2) / (ncol(mat) - 1L)
E_cond_var <- cond_var_acc / n_sub

mu_grid_vec   <- fill_grid(rowMeans(mu_grid_draws))
f_warm_vec    <- fill_grid(rowMeans(f_warm_draws))
spde_cold_vec <- fill_grid(rowMeans(a_cold_draws))
spde_warm_vec <- fill_grid(rowMeans(a_warm_draws))
sd_f_warm_vec <- fill_grid(sqrt(E_cond_var + row_var(f_warm_draws)))
sd_a_warm_vec <- fill_grid(sqrt(E_cond_var + row_var(a_warm_draws)))

mean_s2_e <- mean(theta_post[, "sigma2_e"])
cat(sprintf("\n=== INCERTEZZA A POSTERIORI CONTINUA (ANNO %d, GRIGLIA %dx%d) ===\n",
            years[idx_warm], n_grid, n_grid))
cat(sprintf("SD campo latente f(s)   -> media: %.3f °C | 98° pct: %.3f °C | max bordo: %.3f °C\n",
            mean(sd_f_warm_vec, na.rm = TRUE),
            quantile(sd_f_warm_vec, 0.98, na.rm = TRUE),
            max(sd_f_warm_vec, na.rm = TRUE)))
cat(sprintf("SD anomalia f(s) - m(s) -> media: %.3f °C | range: [%.3f, %.3f] °C\n",
            mean(sd_a_warm_vec, na.rm = TRUE),
            min(sd_a_warm_vec, na.rm = TRUE),
            max(sd_a_warm_vec, na.rm = TRUE)))
cat(sprintf("SD predittiva Y*(s)     -> media: %.3f °C (include sigma_eps = %.3f °C)\n",
            mean(sqrt(sd_f_warm_vec^2 + mean_s2_e), na.rm = TRUE), sqrt(mean_s2_e)))

# ------------------------------------------------------------------------------
# 8. VISUALIZZAZIONI PULITE IN GGPLOT2 (ZERO WARNING, ASSI DINAMICI)
# ------------------------------------------------------------------------------
# Interpolazione P1 lineare sui soli dati osservati per confronto smooth
mesh_obs <- fmesher::fm_rcdt_2d_inla(loc = real_data$coords)
A_obs    <- fmesher::fm_basis(mesh_obs, loc = grid_unit)
in_hull  <- Matrix::rowSums(A_obs) > 0.999

interp_obs_to_grid <- function(vals_obs) {
  v_nodes <- numeric(mesh_obs$n)
  v_nodes[mesh_obs$idx$loc] <- vals_obs
  out <- as.numeric(A_obs %*% v_nodes)
  out[!in_hull] <- NA_real_
  out
}

obs_cold_vec  <- interp_obs_to_grid(anom_emp_obs[, idx_cold])
obs_warm_vec  <- interp_obs_to_grid(anom_emp_obs[, idx_warm])
diff_cold_vec <- obs_cold_vec - spde_cold_vec
diff_warm_vec <- obs_warm_vec - spde_warm_vec

theme_map <- function() {
  theme_minimal(base_size = 10) +
    theme(
      panel.grid       = element_blank(),
      panel.border     = element_rect(color = "grey30", fill = NA, linewidth = 0.5),
      panel.spacing    = unit(0.9, "lines"),
      strip.text       = element_text(face = "bold", size = 9.5, margin = margin(b = 4)),
      plot.title       = element_text(face = "bold", hjust = 0.5, size = 10.5),
      plot.subtitle    = element_text(hjust = 0.5, size = 8.5, color = "grey35"),
      axis.text        = element_text(size = 8.5),
      legend.position  = "bottom",
      legend.title     = element_text(size = 8.5, face = "bold"),
      legend.text      = element_text(size = 8),
      legend.key.width = unit(1.0, "cm"),
      legend.key.height = unit(0.28, "cm"),
      legend.margin    = margin(t = 2),
      plot.margin      = margin(6, 6, 4, 6)
    )
}

# FIGURA A: Fit del campo medio (Scatter osservato vs stimato + Superficie continua m(s))
df_scatter <- data.frame(obs = y_emp_mean, hat = mu_hat_obs)
p_scatter <- ggplot(df_scatter, aes(x = obs, y = hat)) +
  geom_point(color = "royalblue", alpha = 0.7, size = 1.6) +
  geom_abline(slope = 1, intercept = 0, color = "firebrick", linetype = "dashed", linewidth = 0.8) +
  coord_equal() +
  labs(
    title    = "Fit campo medio spaziale",
    subtitle = sprintf("RMSE = %.3f °C (1983–2007)", sqrt(mean((y_emp_mean - mu_hat_obs)^2))),
    x        = "Media osservata (°C)",
    y        = "Stima a posteriori Phi * mu (°C)"
  ) +
  theme_map()

df_mu <- data.frame(base_grid, val = mu_grid_vec)
p_mu <- ggplot(df_mu, aes(x = x, y = y)) +
  geom_raster(aes(fill = val), interpolate = TRUE) +
  geom_contour(aes(z = val), color = "grey15", alpha = 0.45, linewidth = 0.28, bins = 8, na.rm = TRUE) +
  scale_fill_gradientn(colors = hcl.colors(41, "YlOrRd", rev = TRUE),
                       na.value = "grey95", name = "Temp. media (°C)") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3), expand = c(0, 0)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5), expand = c(0, 0)) +
  coord_fixed() +
  labs(
    title    = "Superficie media stimata m(s)",
    subtitle = sprintf("Dominio %.0f km su Roma", scale_km),
    x = "x (km da Roma)", y = "y (km da Roma)"
  ) +
  theme_map()

grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = 1, ncol = 2)))
print(p_scatter, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_mu,      vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))

# Helper per confronto anno singolo: [Osservata P1 | Stimata SPDE] + [Diff. amplificata]
plot_year_comparison <- function(year_label, obs_vec, spde_vec, diff_vec) {
  lim_yr  <- ceiling(max(abs(c(obs_vec, spde_vec)), na.rm = TRUE) * 10) / 10
  lim_df  <- max(0.2, ceiling(quantile(abs(diff_vec), 0.99, na.rm = TRUE) * 10) / 10)
  rmse_yr <- sqrt(mean(diff_vec^2, na.rm = TRUE))

  df_pair <- rbind(
    data.frame(base_grid, val = obs_vec,  panel = sprintf("Osservata %d", year_label)),
    data.frame(base_grid, val = spde_vec, panel = sprintf("Stimata SPDE %d", year_label))
  )
  df_pair$panel <- factor(df_pair$panel, levels = unique(df_pair$panel))

  df_diff <- data.frame(
    base_grid,
    val   = diff_vec,
    panel = sprintf("Diff. Oss. − SPDE (RMS %.2f °C)", rmse_yr)
  )

  p_pair <- ggplot(df_pair, aes(x = x, y = y)) +
    geom_raster(aes(fill = val), interpolate = TRUE) +
    geom_contour(aes(z = val), color = "grey15", alpha = 0.45, linewidth = 0.28, bins = 8, na.rm = TRUE) +
    facet_wrap(~ panel, nrow = 1) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0, limits = c(-lim_yr, lim_yr), oob = scales::squish,
      na.value = "grey95", name = "Anomalia (°C)"
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3), expand = c(0, 0)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5), expand = c(0, 0)) +
    coord_fixed() +
    labs(x = "x (km da Roma)", y = "y (km da Roma)") +
    theme_map()

  p_diff <- ggplot(df_diff, aes(x = x, y = y)) +
    geom_raster(aes(fill = val), interpolate = TRUE) +
    geom_contour(aes(z = val), color = "grey20", alpha = 0.35, linewidth = 0.25, bins = 6, na.rm = TRUE) +
    facet_wrap(~ panel, nrow = 1) +
    scale_fill_gradient2(
      low = "#5E3C99", mid = "#F7F7F7", high = "#E66101",
      midpoint = 0, limits = c(-lim_df, lim_df), oob = scales::squish,
      na.value = "grey95", name = sprintf("Δ °C (±%.1f)", lim_df)
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3), expand = c(0, 0)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5), expand = c(0, 0)) +
    coord_fixed() +
    labs(x = "x (km da Roma)", y = NULL) +
    theme_map()

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(
    nrow = 1, ncol = 2, widths = grid::unit(c(2.05, 1.0), "null")
  )))
  print(p_pair, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_diff, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
}

# FIGURA B & C: Confronto smooth anno più freddo e anno più caldo
plot_year_comparison(years[idx_cold], obs_cold_vec, spde_cold_vec, diff_cold_vec)
plot_year_comparison(years[idx_warm], obs_warm_vec, spde_warm_vec, diff_warm_vec)

# FIGURA D: Mappe di incertezza a posteriori continua (SD)
df_sd  <- data.frame(base_grid, sd_f = sd_f_warm_vec, sd_a = sd_a_warm_vec)
obs_df <- data.frame(x = coords_km[, 1], y = coords_km[, 2])

cap_f <- as.numeric(quantile(df_sd$sd_f, 0.98, na.rm = TRUE))
lim_f <- c(min(df_sd$sd_f, na.rm = TRUE), cap_f)
lim_a <- range(df_sd$sd_a, na.rm = TRUE)

p_sd_f <- ggplot(df_sd, aes(x = x, y = y)) +
  geom_raster(aes(fill = sd_f), interpolate = TRUE) +
  geom_contour(aes(z = sd_f), color = "white", alpha = 0.4, linewidth = 0.22,
               breaks = pretty(lim_f, n = 5), na.rm = TRUE) +
  geom_point(data = obs_df, aes(x = x, y = y), size = 0.3, alpha = 0.25, color = "black") +
  scale_fill_viridis_c(
    option = "magma", direction = -1,
    limits = lim_f, oob = scales::squish, na.value = "grey95", name = "SD (°C)"
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3), expand = c(0, 0)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5), expand = c(0, 0)) +
  coord_fixed() +
  labs(
    title    = sprintf("SD campo totale f_%d(s)", years[idx_warm]),
    subtitle = sprintf("Media %.3f °C (scala 98° pct: %.2f °C)", mean(df_sd$sd_f, na.rm = TRUE), cap_f),
    x = "x (km da Roma)", y = "y (km da Roma)"
  ) +
  theme_map()

p_sd_a <- ggplot(df_sd, aes(x = x, y = y)) +
  geom_raster(aes(fill = sd_a), interpolate = TRUE) +
  geom_contour(aes(z = sd_a), color = "white", alpha = 0.4, linewidth = 0.22, bins = 6, na.rm = TRUE) +
  geom_point(data = obs_df, aes(x = x, y = y), size = 0.3, alpha = 0.25, color = "black") +
  scale_fill_viridis_c(
    option = "magma", direction = -1,
    limits = lim_a, na.value = "grey95", name = "SD (°C)"
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3), expand = c(0, 0)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5), expand = c(0, 0)) +
  coord_fixed() +
  labs(
    title    = sprintf("SD anomalia %d (f − m)", years[idx_warm]),
    subtitle = sprintf("Media %.3f °C (range %.2f–%.2f °C)", mean(df_sd$sd_a, na.rm = TRUE), lim_a[1], lim_a[2]),
    x = "x (km da Roma)", y = NULL
  ) +
  theme_map()

grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = 1, ncol = 2)))
print(p_sd_f, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_sd_a, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
