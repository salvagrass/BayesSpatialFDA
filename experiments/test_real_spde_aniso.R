# ==============================================================================
# REAL-DATA PLANAR TEST - ANISOTROPIC SPDE PC-GIBBS SAMPLER (C++ RCPP BACKEND)
# Model: Y_i = Phi %*% beta_i + eps_i,   beta_i | mu, theta ~ N_V(mu, Q(theta)^-1)
# ==============================================================================

devtools::load_all(compile = FALSE)
library(Matrix)
library(fmesher)
library(coda)
library(ggplot2)
library(grid)

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

Y_mat    <- as.matrix(real_data$y_mat)
Phi      <- as(real_data$Phi, "dgCMatrix")
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
    "Domain boundary",
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

Gx    <- as(spde_matrices$K11, "dgCMatrix")
Gy    <- as(spde_matrices$K22, "dgCMatrix")
Gxy   <- as(spde_matrices$K12, "dgCMatrix")
C_vec <- as.numeric(Matrix::diag(spde_matrices$C))

# ------------------------------------------------------------------------------
# 3. ESECUZIONE GIBBS SAMPLER C++ (MULTI-CHAIN)
# ------------------------------------------------------------------------------
num_iteration <- 5000L
num_burn_in   <- 1500L

cat("\n[R] Avvio C++ Gibbs Engine sui dati reali...\n")
t_start <- proc.time()

fit <- .run_gibbs_sampler_cpp(
  Phi         = Phi,
  Y           = Y_mat,
  Gx          = Gx,
  Gy          = Gy,
  Gxy         = Gxy,
  C           = C_vec,
  num_iters   = num_iteration,
  num_burn_in = num_burn_in
)

elapsed    <- (proc.time() - t_start)[["elapsed"]]
num_chains <- length(fit)
cat(sprintf(
  "\n[R] C++ Sampler completato in %.2f s | Catene: %d | Iterazioni/catena: %d (burn-in: %d)\n",
  elapsed, num_chains, num_iteration, num_burn_in
))

# ------------------------------------------------------------------------------
# 4. DIAGNOSTICA MCMC (CODA, ESS, R-HAT) E RIEPILOGO FISICO A POSTERIORI
# ------------------------------------------------------------------------------
chain_draws_list <- lapply(fit, function(chain) {
  th_mat <- do.call(rbind, chain$theta)
  colnames(th_mat) <- c("kappa", "sigma_beta", "v1", "v2")
  s2_vec <- as.numeric(unlist(chain$sigma_epsilon))
  cbind(th_mat, sigma2_e = s2_vec)
})

mcmc_post <- coda::mcmc.list(lapply(chain_draws_list, function(mat) {
  coda::mcmc(mat, start = num_burn_in + 1L)
}))

theta_post <- do.call(rbind, chain_draws_list)

# Traceplot multi-chain in scala fisica interpretabile (post-burn-in)
chain_cols <- adjustcolor(c("grey25", "steelblue", "firebrick", "forestgreen")[seq_len(num_chains)], 0.75)

plot_multi_trace <- function(vals_list, main, ylab) {
  all_v <- unlist(vals_list)
  yl    <- range(all_v, na.rm = TRUE)
  n_s   <- length(vals_list[[1]])
  x_it  <- seq.int(num_burn_in + 1L, num_burn_in + n_s)
  plot(x_it, vals_list[[1]], type = "l", col = chain_cols[1],
       xlab = "Iterazione (post burn-in)", ylab = ylab, main = main, ylim = yl)
  if (length(vals_list) > 1L) {
    for (ch in 2:length(vals_list)) {
      lines(x_it, vals_list[[ch]], col = chain_cols[ch])
    }
  }
}

par(mfrow = c(2, 3), mar = c(4, 4, 2.2, 1))
plot_multi_trace(lapply(chain_draws_list, function(m) sqrt(m[, "sigma2_e"])),
                 "sigma_eps (°C)", "°C")
plot_multi_trace(lapply(chain_draws_list, function(m) (sqrt(8) / m[, "kappa"]) * scale_km),
                 "Range rho (km)", "km")
plot_multi_trace(lapply(chain_draws_list, function(m) m[, "sigma_beta"]),
                 "sigma_beta (°C)", "°C")
plot_multi_trace(lapply(chain_draws_list, function(m) m[, "v1"]),
                 "Anisotropia v1", "v1")
plot_multi_trace(lapply(chain_draws_list, function(m) m[, "v2"]),
                 "Anisotropia v2", "v2")
plot_multi_trace(lapply(chain_draws_list, function(m) sqrt(m[, "v1"]^2 + m[, "v2"]^2)),
                 "Norma ||v||", "||v||")
par(mfrow = c(1, 1))

# Riepilogo coda + ESS + Gelman-Rubin
print(summary(mcmc_post))
cat("\nEffective Sample Size (ESS totale):\n")
print(round(coda::effectiveSize(mcmc_post), 1))
if (num_chains > 1L) {
  cat("\nGelman-Rubin diagnostic (R-hat):\n")
  print(coda::gelman.diag(mcmc_post, autoburnin = FALSE))
}

# Grandezze fisiche derivate
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
# 5. CAMPO MEDIO m(s) E ANOMALIE (MULTI-CHAIN POST-BURN-IN)
# ------------------------------------------------------------------------------
mu_samples_mat <- do.call(rbind, lapply(fit, function(chain) do.call(rbind, chain$mu)))
mu_post_draws  <- t(mu_samples_mat) # V x n_post
mu_post_mean   <- rowMeans(mu_post_draws)
mu_hat_obs     <- as.numeric(Phi %*% mu_post_mean)
y_emp_mean     <- rowMeans(Y_mat)
anom_emp_obs   <- sweep(Y_mat, 1, y_emp_mean, "-")

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
# 6. PREDIZIONE CONTINUA (80x80) E INCERTEZZA A POSTERIORI (VARIANZA TOTALE)
# ------------------------------------------------------------------------------
loc_label <- "Roma" # Aggiorna quando cambi località in Sezione 1

n_grid <- 80L
x_seq_unit <- seq(min(real_data$boundary[, 1]), max(real_data$boundary[, 1]), length.out = n_grid)
y_seq_unit <- seq(min(real_data$boundary[, 2]), max(real_data$boundary[, 2]), length.out = n_grid)
grid_unit  <- as.matrix(expand.grid(x = x_seq_unit, y = y_seq_unit))
base_grid  <- expand.grid(x = x_seq_unit * scale_km, y = y_seq_unit * scale_km)

# Triangolazione di Delaunay pura sulle sole N osservazioni (0 nodi extra a 0)
mesh_obs <- if (exists("fm_delaunay_2d", envir = asNamespace("fmesher"))) {
  fmesher::fm_delaunay_2d(real_data$coords)
} else {
  fmesher::fm_rcdt_2d_inla(loc = real_data$coords)
}
A_obs     <- fmesher::fm_basis(mesh_obs, loc = grid_unit)
extra_obs <- setdiff(seq_len(mesh_obs$n), mesh_obs$idx$loc)
in_hull   <- Matrix::rowSums(A_obs) > 0.999
if (length(extra_obs) > 0L) {
  in_hull <- in_hull & (Matrix::rowSums(A_obs[, extra_obs, drop = FALSE]) < 1e-8)
}

# Valutazione FEM continua escludendo anche i pixel agganciati ai nodi vuoti del buffer
Phi_grid_all <- as.matrix(fmesher::fm_basis(mesh, loc = grid_unit))
in_mesh      <- in_hull
if (length(empty_idx) > 0L) {
  in_mesh <- in_mesh & (rowSums(Phi_grid_all[, empty_idx, drop = FALSE]) < 1e-6)
}

Phi_grid_in <- Phi_grid_all[in_mesh, , drop = FALSE]
n_in        <- nrow(Phi_grid_in)

# Allineamento esatto catena per catena tra mu (thinnato) e theta/sigma2_e
aligned_mu_list <- vector("list", num_chains)
aligned_th_list <- vector("list", num_chains)
for (ch in seq_len(num_chains)) {
  mu_ch  <- do.call(rbind, fit[[ch]]$mu)
  th_ch  <- chain_draws_list[[ch]]
  n_m    <- nrow(mu_ch)
  n_t    <- nrow(th_ch)
  th_idx <- round(seq(n_t / n_m, n_t, length.out = n_m))
  aligned_mu_list[[ch]] <- t(mu_ch)
  aligned_th_list[[ch]] <- th_ch[th_idx, , drop = FALSE]
}
mu_aligned <- do.call(cbind, aligned_mu_list) # V x n_sub
th_aligned <- do.call(rbind, aligned_th_list) # n_sub x 5
n_sub      <- nrow(th_aligned)

PtP <- Matrix::crossprod(Phi)
PtY <- as.matrix(Matrix::crossprod(Phi, Y_mat))

cond_var_acc  <- numeric(n_in)
mu_grid_draws <- matrix(0, nrow = n_in, ncol = n_sub)
a_cold_draws  <- matrix(0, nrow = n_in, ncol = n_sub)
a_warm_draws  <- matrix(0, nrow = n_in, ncol = n_sub)
f_warm_draws  <- matrix(0, nrow = n_in, ncol = n_sub)
beta_all_acc  <- matrix(0, nrow = V, ncol = M)

for (j in seq_len(n_sub)) {
  mu_s <- mu_aligned[, j]

  Q_th <- build_precision_matrix(
    spde_matrices = spde_matrices,
    kappa         = th_aligned[j, "kappa"],
    sigma_beta    = th_aligned[j, "sigma_beta"],
    v1            = th_aligned[j, "v1"],
    v2            = th_aligned[j, "v2"]
  )
  inv_s2 <- 1 / th_aligned[j, "sigma2_e"]

  Q_post     <- Matrix::forceSymmetric(inv_s2 * PtP + Q_th)
  Sigma_post <- as.matrix(Matrix::solve(Q_post))

  # Varianza condizionale di phi(s)^T beta_i sui punti interni
  cond_var_acc <- cond_var_acc + rowSums((Phi_grid_in %*% Sigma_post) * Phi_grid_in)

  # Medie condizionali di tutte le anomalie (V x M)
  rhs_all      <- inv_s2 * sweep(PtY, 1, as.numeric(PtP %*% mu_s), "-")
  m_anom       <- Sigma_post %*% rhs_all
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
cat(sprintf("In-sample RMSE campo totale -> tutti gli anni: %.3f °C | anno 2007: %.3f °C (n_sub = %d draw allineati)\n",
            rmse_all, rmse_2007, n_sub))

# Ricomposizione vettori completi (NA fuori dal dominio utile)
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
# 7. VISUALIZZAZIONI PULITE IN GGPLOT2 (ZERO WARNING, ASSI DINAMICI)
# ------------------------------------------------------------------------------
interp_obs_to_grid <- function(vals_obs) {
  v_nodes <- numeric(mesh_obs$n)
  v_nodes[mesh_obs$idx$loc] <- vals_obs
  out <- as.numeric(A_obs %*% v_nodes)
  out[!in_mesh] <- NA_real_
  out
}

obs_cold_vec  <- interp_obs_to_grid(anom_emp_obs[, idx_cold])
obs_warm_vec  <- interp_obs_to_grid(anom_emp_obs[, idx_warm])
diff_cold_vec <- obs_cold_vec - spde_cold_vec
diff_warm_vec <- obs_warm_vec - spde_warm_vec

theme_map <- function() {
  theme_minimal(base_size = 10) +
    theme(
      panel.grid        = element_blank(),
      panel.border      = element_rect(color = "grey30", fill = NA, linewidth = 0.5),
      panel.spacing     = unit(0.9, "lines"),
      strip.text        = element_text(face = "bold", size = 9.5, margin = margin(b = 4)),
      plot.title        = element_text(face = "bold", hjust = 0.5, size = 10.5),
      plot.subtitle     = element_text(hjust = 0.5, size = 8.5, color = "grey35"),
      axis.text         = element_text(size = 8.5),
      legend.position   = "bottom",
      legend.title      = element_text(size = 8.5, face = "bold"),
      legend.text       = element_text(size = 8),
      legend.key.width  = unit(1.0, "cm"),
      legend.key.height = unit(0.28, "cm"),
      legend.margin     = margin(t = 2),
      plot.margin       = margin(6, 6, 4, 6)
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
    subtitle = sprintf("Dominio %.0f km su %s", scale_km, loc_label),
    x = sprintf("x (km da %s)", loc_label), y = sprintf("y (km da %s)", loc_label)
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
    labs(x = sprintf("x (km da %s)", loc_label), y = sprintf("y (km da %s)", loc_label)) +
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
    labs(x = sprintf("x (km da %s)", loc_label), y = NULL) +
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
    x = sprintf("x (km da %s)", loc_label), y = sprintf("y (km da %s)", loc_label)
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
    x = sprintf("x (km da %s)", loc_label), y = NULL
  ) +
  theme_map()

grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = 1, ncol = 2)))
print(p_sd_f, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_sd_a, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
