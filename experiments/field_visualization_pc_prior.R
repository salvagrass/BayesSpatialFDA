# Prior predictive check for the PC prior on theta = (v1, v2, kappa, sigma_beta)
# of the ANISOTROPIC SPDE, as realized fields.

# The PC prior is not set through its rates (lambda_v, lambda_theta,
# lambda_sigma_beta) but through three tail statements:
#   P(exp(||v||) > a0)   = alpha_v      (a0: surprisingly strong anisotropy)
#   P(rho < rho0)        = alpha_rho    (rho0: surprisingly short range)
#   P(sigma_beta > sigma0) = alpha_sigma  (sigma0: surprisingly large sd)
# with rho = sqrt(8) / kappa the correlation range.
# This script takes a choice of (a0, rho0, sigma0) and gives plots
# of fields generated according to those choices.
#
# The map Psi pushes a standard Gaussian Y ~ N(0, I_4) forward to the PC prior
# exactly, so drawing from the prior is: draw Y, apply psi_reparam_vec().
#
# Four panels: a baseline at the prior medians of (||v||, kappa, sigma_beta),
# and three panels that each move ONE of those to a rare value (its TAIL_Q
# prior quantile) while the others stay at the baseline. The changed parameter
# is highlighted in the panel title. The direction alpha = arg(v)/2 is the
# same in all four.
# Every panel uses the SAME random seed (common random numbers), so differences
# between panels come from the parameters, not from fresh noise. Each panel has
# its own colour scale and its own legend: the same colour means different
# values in different panels, so read the amplitude from the legend. On a
# bounded mesh a range much longer than the domain shows up as a nearly flat
# field with a large offset.

library(Matrix)
library(fmesher)
devtools::load_all(quiet = TRUE)

# ------------------------------------------------------------------------------
#   * rho0 (range) is a LENGTH, so it only means something relative to the
#     domain.
#   * a0 (anisotropy ratio exp(||v||)) is DIMENSIONLESS.
#   * sigma0 (marginal sd) is in the UNITS OF THE RESPONSE.
# The domain is therefore defined first, and the range statement is tied to it.
MESH_SIDE <- 10 # domain is [0, MESH_SIDE] x [0, MESH_SIDE]
MESH_N <- 90 # reference mesh built from a MESH_N x MESH_N grid of points

# ------------------------------------------------------------------------------
a0          <- 10    # anisotropy ratio exp(||v||) considered surprisingly large
alpha_v     <- 0.05  # P(exp(||v||) > a0) = alpha_v

# Range statement, as a fraction of the domain side.
# rho0 = MESH_SIDE/10 follows the simulation anchor of Llamazares-Elias et al.
# Shrink rho0 for a finer baseline (needs a denser mesh), grow it for a smoother one.
rho0        <- 0.10 * MESH_SIDE  # correlation range considered surprisingly short
alpha_rho   <- 0.05  # P(rho < rho0) = alpha_rho

sigma0      <- 2.0   # marginal sd considered surprisingly large
alpha_sigma <- 0.05  # P(sigma_beta > sigma0) = alpha_sigma

# ------------------------------------------------------------------------------
SEED <- 2024 # seed of the prior draws and of the common random numbers
TAIL_Q <- 0.99 # prior quantile that defines a "rare" value in the field panels
ALPHA_DEG <- 30 # direction arg(v)/2 of the strongest correlation, all panels
FIELD_RES <- 400 # pixels per side of each field panel
FIELDS_PNG <- NULL # file path: also save the field panels as a slide-ready PNG

# ------------------------------------------------------------------------------
# CALIBRATION
# ------------------------------------------------------------------------------
# pc_calibrate_all() takes the range statement on the kappa scale. Since
# rho = sqrt(8) / kappa, P(rho < rho0) = P(kappa > kappa0) with
# kappa0 = sqrt(8) / rho0.
kappa0 <- sqrt(8) / rho0
lambdas <- pc_calibrate_all(a0, alpha_v, kappa0, alpha_rho, sigma0, alpha_sigma)

cat("\n=== PC-PRIOR RATES ===\n")
cat(sprintf(
  "lambda_v = %.4f | lambda_theta = %.4f | lambda_sigma_beta = %.4f\n",
  lambdas$lambda_v,
  lambdas$lambda_theta,
  lambdas$lambda_sigma_beta
))

# ------------------------------------------------------------------------------
# REALIZED FIELDS
# ------------------------------------------------------------------------------
# Baseline: the prior medians of kappa, ||v|| and sigma_beta. Each of the other
# three panels moves ONE of them to its TAIL_Q prior quantile (short range,
# strong anisotropy, large sd) and keeps the other two at the baseline. Medians
# and quantiles are those of the marginal priors. The
# direction alpha = arg(v)/2 is uniform under the prior, so no value is more
# typical than another: it is fixed at ALPHA_DEG in every panel.
x <- seq(0, MESH_SIDE, length.out = MESH_N)
y <- seq(0, MESH_SIDE, length.out = MESH_N)
locs <- as.matrix(expand.grid(x = x, y = y))
mesh <- fm_mesh_2d(loc = locs)

baseline <- c(
  kappa = median(kappa),
  norm_v = median(norm_v),
  sigma_beta = median(sigma_beta)
)
rare <- c(
  kappa = unname(quantile(kappa, TAIL_Q)),
  norm_v = unname(quantile(norm_v, TAIL_Q)),
  sigma_beta = unname(quantile(sigma_beta, TAIL_Q))
)

field_tbl <- data.frame(
  panel = c(
    "Baseline (prior medians)",
    "Rare: short range",
    "Rare: strong anisotropy",
    "Rare: large amplitude"
  ),
  changed = c(NA, "kappa", "norm_v", "sigma_beta"),
  kappa = baseline[["kappa"]],
  norm_v = baseline[["norm_v"]],
  sigma_beta = baseline[["sigma_beta"]]
)
for (k in which(!is.na(field_tbl$changed))) {
  field_tbl[k, field_tbl$changed[k]] <- rare[[field_tbl$changed[k]]]
}
alpha <- ALPHA_DEG * pi / 180
field_tbl$v1 <- field_tbl$norm_v * cos(2 * alpha)
field_tbl$v2 <- field_tbl$norm_v * sin(2 * alpha)

fields <- lapply(seq_len(nrow(field_tbl)), function(k) {
  p <- field_tbl[k, ]
  sim <- generate_toy_data_pc_prior(
    mesh = mesh,
    num_year = 1,
    mu_prior_mean = 0.0,
    Q_mu = 1.0,
    kappa = p$kappa,
    v1 = p$v1,
    v2 = p$v2,
    sigma_beta = p$sigma_beta,
    sigma_e = 1e-6,
    seed = SEED # same seed in every panel: common random numbers
  )
  sim$beta_mat[, 1] - sim$mu_true # latent anomaly
})

# Panel title: the scenario, then the four parameters, with the one that
# differs from the baseline in colour.
draw_field_title <- function(p) {
  terms <- list(
    kappa = bquote(kappa == .(sprintf("%.2f", p$kappa))),
    norm_v = bquote(
      group("|", group("|", bold(v), "|"), "|") == .(sprintf("%.2f", p$norm_v))
    ),
    sigma_beta = bquote(sigma[beta] == .(sprintf("%.2f", p$sigma_beta))),
    alpha = bquote(alpha == .(ALPHA_DEG) * degree)
  )
  # The parameter line with only the terms flagged in `show` visible; the
  # others keep their space, so two passes give a two-colour line.
  param_line <- function(show) {
    shown <- Map(
      function(term, visible) if (visible) term else call("phantom", term),
      terms,
      show
    )
    as.expression(Reduce(function(a, b) bquote(.(a) * "     " * .(b)), shown))
  }
  is_changed <- names(terms) %in% p$changed

  x_mid <- mean(range(mesh$loc[, 1]))
  y_top <- max(mesh$loc[, 2])
  y_par <- y_top + 1.3 * strheight("M")
  y_name <- y_top + 3.3 * strheight("M")
  text(x_mid, y_name, p$panel, font = 2, cex = 1.15, xpd = NA)
  text(x_mid, y_par, param_line(!is_changed), col = "grey25", xpd = NA)
  if (any(is_changed)) {
    text(x_mid, y_par, param_line(is_changed), col = "#C8102E", xpd = NA)
  }
}

draw_field_panels <- function() {
  old_par <- par(mfrow = c(2, 2), mar = c(2, 2.4, 3.6, 0.4))
  on.exit(par(old_par))
  for (k in seq_along(fields)) {
    plot_spde_field(mesh, fields[[k]], main_title = "", resolution = FIELD_RES)
    draw_field_title(field_tbl[k, ])
  }
}

draw_field_panels()
if (!is.null(FIELDS_PNG)) {
  png(FIELDS_PNG, width = 6.6, height = 5.8, units = "in", res = 400)
  draw_field_panels()
  dev.off()
}

cat("\n=== FIELD CONFIGURATIONS (one row per panel) ===\n")
print(
  data.frame(
    panel = field_tbl$panel,
    kappa = round(field_tbl$kappa, 3),
    rho = round(sqrt(8) / field_tbl$kappa, 2),
    norm_v = round(field_tbl$norm_v, 3),
    aniso_ratio = round(exp(field_tbl$norm_v), 2),
    sigma_beta = round(field_tbl$sigma_beta, 3),
    alpha_deg = ALPHA_DEG,
    field_sd = round(vapply(fields, sd, numeric(1)), 3)
  ),
  row.names = FALSE
)
