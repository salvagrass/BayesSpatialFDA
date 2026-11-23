# Diagnostic plots and accuracy metrics for the PC-parameterized Gibbs sampler
# (experiments/pc_gibbs/run_sampler.R, branch feature/sampler-cpp).
#
# NOTE: this file is written for the output of the pc_gibbs sampler
# (modeling_choices_v2 parameterization: kappa, sigma_beta, v1, v2). It does NOT
# work with the older sampler in experiments/gibbs_sampler/, which has its own
# diagnostics.R with functions of the same names: source only one of the two
# files in a given R session.
#
# Requirements: packages coda, ggplot2, patchwork, fmesher, Matrix.
# Typical use, at the end of experiments/pc_gibbs_prototype.R:
#   source("experiments/pc_gibbs/diagnostics.R")
#   diag_results <- run_all_diagnostics(
#     res = fit, mesh = mesh, sim_data = sim_data, truth = truth,
#     num_burn_in = num_burn_in, elapsed = elapsed, year = 1L, vertex = 1L
#   )
#
# Conventions of run_gibbs_sampler() (pc_gibbs):
#   theta_samples  4 x num_iteration (kappa, sigma_beta, v1, v2), full chain
#   sigma2_samples num_iteration, full chain
#   beta_samples   V x M x num_saved, post burn-in, CENTERED: beta_i = mu + beta_tilde_i
#   mu_samples     V x num_saved, post burn-in
#   acceptance_rate
# The anomaly is therefore beta_tilde_i = beta_i - mu, draw by draw.
#
# `sim_data` is the output of a toy data generator of the sampler-cpp family
# (generate_toy_data, generate_toy_data_matern, generate_toy_data_nonstationary):
# it must contain y_mat, mu_true, beta_mat, Phi and locs.
# `truth` is the list used to generate the data: kappa, sigma_beta, v1, v2, sigma_e.
# For a misspecified generator (e.g. generate_toy_data_matern with nu != 1) the
# true kappa is only a "same range" reference, not a value the model can recover.

#' Posterior Summaries of the Anomaly Field beta_tilde = beta - mu
#'
#' @param res Output list returned by `run_gibbs_sampler()`
#' @return List with `mean` and `sd` (V x M matrices)
summarize_anomaly <- function(res) {
  dims <- dim(res$beta_samples)
  V <- dims[1]
  M <- dims[2]
  anomaly_mean <- matrix(NA_real_, V, M)
  anomaly_sd <- matrix(NA_real_, V, M)
  for (t in seq_len(M)) {
    draws_t <- res$beta_samples[, t, ] - res$mu_samples
    anomaly_mean[, t] <- rowMeans(draws_t)
    anomaly_sd[, t] <- apply(draws_t, 1, sd)
  }
  list(mean = anomaly_mean, sd = anomaly_sd)
}

#' True Anomaly Field at the Mesh Vertices
#'
#' @param sim_data Output list of the toy data generator
#' @return V x M matrix beta - mu
get_true_anomaly <- function(sim_data) {
  sim_data$beta_mat - sim_data$mu_true
}

#' Common ggplot2 Theme for the Diagnostic Plots
#'
#' @param base_size Base font size
diag_theme <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey92", colour = NA),
      panel.grid.minor = ggplot2::element_blank()
    )
}

#' Plot Full MCMC Traceplots with Burn-in and Ground Truth (ggplot2)
#'
#' Hyperparameters are stored for the full chain; beta_tilde[vertex, 1] only
#' after burn-in, so its trace starts at `num_burn_in + 1`. The y-axis of each
#' panel is zoomed on the draws after half of the burn-in (plus the true value):
#' the first iterations can be orders of magnitude away (e.g. sigma^2_e starts
#' around 100) and would otherwise flatten the chain. Earlier values outside
#' that window are drawn at the edge of the panel.
#'
#' @param res Output list returned by `run_gibbs_sampler()`
#' @param sim_data Output list of the toy data generator
#' @param truth List of true values (kappa, sigma_beta, v1, v2, sigma_e)
#' @param num_burn_in Number of burn-in iterations
#' @param vertex Vertex index used for the beta_tilde trace
#' @return Invisibly returns the ggplot object (it is also printed)
plot_mcmc_traces <- function(res, sim_data, truth, num_burn_in, vertex = 1L) {
  th <- res$theta_samples
  num_iteration <- ncol(th)
  true_anomaly <- get_true_anomaly(sim_data)
  beta_label <- sprintf("beta_tilde[%d,1]", vertex)

  traces <- list(
    list(name = "sigma^2_e", iter = seq_len(num_iteration),
         value = res$sigma2_samples, true = truth$sigma_e^2),
    list(name = "kappa", iter = seq_len(num_iteration),
         value = th["kappa", ], true = truth$kappa),
    list(name = "sigma_beta", iter = seq_len(num_iteration),
         value = th["sigma_beta", ], true = truth$sigma_beta),
    list(name = "v1", iter = seq_len(num_iteration),
         value = th["v1", ], true = truth$v1),
    list(name = "v2", iter = seq_len(num_iteration),
         value = th["v2", ], true = truth$v2),
    list(name = beta_label, iter = (num_burn_in + 1L):num_iteration,
         value = res$beta_samples[vertex, 1, ] - res$mu_samples[vertex, ],
         true = true_anomaly[vertex, 1])
  )
  param_levels <- vapply(traces, `[[`, character(1), "name")

  zoom_from <- floor(num_burn_in / 2)
  trace_df <- do.call(rbind, lapply(traces, function(tr) {
    keep <- tr$iter > zoom_from
    lim <- range(c(tr$value[keep], tr$true))
    pad <- 0.05 * diff(lim) + 1e-12
    data.frame(
      parameter = tr$name,
      iteration = tr$iter,
      value = pmin(pmax(tr$value, lim[1] - pad), lim[2] + pad)
    )
  }))
  trace_df$parameter <- factor(trace_df$parameter, levels = param_levels)

  truth_df <- data.frame(
    parameter = factor(param_levels, levels = param_levels),
    true = vapply(traces, `[[`, numeric(1), "true")
  )

  p <- ggplot2::ggplot(trace_df, ggplot2::aes(x = iteration, y = value)) +
    ggplot2::geom_line(linewidth = 0.25, colour = "grey20") +
    ggplot2::geom_vline(xintercept = num_burn_in, linetype = "dashed", colour = "steelblue") +
    ggplot2::geom_hline(
      data = truth_df,
      ggplot2::aes(yintercept = true),
      colour = "firebrick",
      linewidth = 0.8
    ) +
    ggplot2::facet_wrap(~parameter, scales = "free_y", ncol = 3) +
    ggplot2::labs(
      title = "MCMC traces",
      subtitle = "Red: true value | dashed blue: end of burn-in | y-axis zoomed from half of the burn-in",
      x = "Iteration",
      y = NULL
    ) +
    diag_theme()

  print(p)
  invisible(p)
}

#' Build Post Burn-In Matrix on the Scale of the Metropolis Step
#'
#' @param res Output list returned by `run_gibbs_sampler()`
#' @param num_burn_in Number of burn-in iterations
#' @param vertex Vertex index used for beta_tilde[vertex, 1]
#' @return Numeric matrix of post burn-in draws
get_post_draws_matrix <- function(res, num_burn_in, vertex = 1L) {
  idx <- (num_burn_in + 1L):ncol(res$theta_samples)
  th <- res$theta_samples

  draws <- cbind(
    "log sigma^2_e" = log(res$sigma2_samples[idx]),
    "log kappa" = log(th["kappa", idx]),
    "log sigma_beta" = log(th["sigma_beta", idx]),
    "v1" = th["v1", idx],
    "v2" = th["v2", idx],
    "beta_tilde" = res$beta_samples[vertex, 1, ] - res$mu_samples[vertex, ]
  )
  colnames(draws)[6] <- sprintf("beta_tilde[%d,1]", vertex)
  draws
}

#' Plot ACF (ggplot2) and Compute Effective Sample Size (ESS)
#'
#' @param res Output list returned by `run_gibbs_sampler()`
#' @param num_burn_in Number of burn-in iterations
#' @param lag_max Maximum lag for ACF plots
#' @param vertex Vertex index used for beta_tilde[vertex, 1]
#' @return Data frame with ESS and total post-burn-in iterations
plot_mcmc_acf_ess <- function(res, num_burn_in, lag_max = 60L, vertex = 1L) {
  post_draws <- get_post_draws_matrix(res, num_burn_in, vertex)
  n_draws <- nrow(post_draws)

  acf_df <- do.call(rbind, lapply(colnames(post_draws), function(nm) {
    a <- stats::acf(post_draws[, nm], lag.max = lag_max, plot = FALSE)
    data.frame(parameter = nm, lag = as.vector(a$lag), acf = as.vector(a$acf))
  }))
  acf_df$parameter <- factor(acf_df$parameter, levels = colnames(post_draws))
  conf_band <- stats::qnorm(0.975) / sqrt(n_draws)

  p <- ggplot2::ggplot(acf_df, ggplot2::aes(x = lag, y = acf)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
    ggplot2::geom_hline(
      yintercept = c(-conf_band, conf_band),
      linetype = "dashed",
      colour = "steelblue"
    ) +
    ggplot2::geom_segment(ggplot2::aes(xend = lag, yend = 0), linewidth = 0.4) +
    ggplot2::facet_wrap(~parameter, ncol = 3) +
    ggplot2::coord_cartesian(ylim = c(min(-0.2, min(acf_df$acf)), 1)) +
    ggplot2::labs(
      title = "Autocorrelation (post burn-in)",
      subtitle = "Dashed blue: 95% band for zero autocorrelation",
      x = "Lag",
      y = "ACF"
    ) +
    diag_theme()
  print(p)

  ess <- coda::effectiveSize(coda::mcmc(post_draws))
  data.frame(
    Parameter = names(ess),
    ESS = round(as.numeric(ess)),
    Out_of = n_draws,
    row.names = NULL
  )
}

#' Posterior Correlation Matrix of the Hyperparameters
#'
#' @param res Output list returned by `run_gibbs_sampler()`
#' @param num_burn_in Number of burn-in iterations
#' @return Rounded correlation matrix of the 5 hyperparameters
compute_hyper_correlations <- function(res, num_burn_in) {
  post_draws <- get_post_draws_matrix(res, num_burn_in)[, 1:5, drop = FALSE]
  round(cor(post_draws), 2)
}

#' Summarize Posterior Estimates vs Ground Truth
#'
#' Also reports the derived physical quantities of modeling_choices_v2:
#' range rho = sqrt(8) / kappa and anisotropy ratio a = exp(||v||).
#'
#' @param res Output list returned by `run_gibbs_sampler()`
#' @param sim_data Output list of the toy data generator
#' @param truth List of true values (kappa, sigma_beta, v1, v2, sigma_e)
#' @param num_burn_in Number of burn-in iterations
#' @param vertex Vertex index used for beta_tilde and mu
#' @return Matrix with True, Mean, 2.5%, Median, and 97.5% quantiles
summarize_hyperparameters <- function(res, sim_data, truth, num_burn_in, vertex = 1L) {
  idx <- (num_burn_in + 1L):ncol(res$theta_samples)
  th <- res$theta_samples
  true_anomaly <- get_true_anomaly(sim_data)

  summ_vec <- function(x, true_val) {
    q <- quantile(x, c(0.025, 0.5, 0.975))
    c(
      True = if (is.null(true_val)) NA_real_ else true_val,
      Mean = mean(x),
      `2.5%` = q[[1]],
      Median = q[[2]],
      `97.5%` = q[[3]]
    )
  }

  kappa <- th["kappa", idx]
  v_norm <- sqrt(th["v1", idx]^2 + th["v2", idx]^2)

  est_tab <- rbind(
    "sigma^2_e" = summ_vec(res$sigma2_samples[idx], truth$sigma_e^2),
    "kappa" = summ_vec(kappa, truth$kappa),
    "sigma_beta" = summ_vec(th["sigma_beta", idx], truth$sigma_beta),
    "v1" = summ_vec(th["v1", idx], truth$v1),
    "v2" = summ_vec(th["v2", idx], truth$v2),
    "range sqrt(8)/kappa" = summ_vec(sqrt(8) / kappa, sqrt(8) / truth$kappa),
    "aniso ratio exp(|v|)" = summ_vec(exp(v_norm), exp(sqrt(truth$v1^2 + truth$v2^2))),
    "beta_tilde[1,1]" = summ_vec(
      res$beta_samples[vertex, 1, ] - res$mu_samples[vertex, ],
      true_anomaly[vertex, 1]
    ),
    "mu" = summ_vec(res$mu_samples[vertex, ], sim_data$mu_true[vertex])
  )
  rownames(est_tab)[8:9] <- c(sprintf("beta_tilde[%d,1]", vertex), sprintf("mu[%d]", vertex))
  round(est_tab, 4)
}

#' Evaluate Field Reconstruction Accuracy vs Naive Baseline
#'
#' The model is evaluated at the mesh vertices and at the observation
#' locations (through Phi). The naive baseline (each year minus the location's
#' average over years) only exists at the observation locations, so the
#' Model-vs-Naive comparison is made there.
#'
#' @param res Output list returned by `run_gibbs_sampler()`
#' @param sim_data Output list of the toy data generator (y_mat, beta_mat, mu_true, Phi)
#' @param anomaly Output of `summarize_anomaly(res)`
#' @param plot_results Logical; whether to plot the model and naive scatter plots
#' @return Data frame comparing Model vs Naive baseline metrics
evaluate_field_accuracy <- function(res, sim_data, anomaly, plot_results = TRUE) {
  Phi <- sim_data$Phi
  Y_mat <- sim_data$y_mat
  true_vert <- get_true_anomaly(sim_data)

  metrics <- function(est, truth) {
    err <- est - truth
    c(
      sqrt(mean(err^2)),
      mean(abs(err)),
      cor(as.vector(est), as.vector(truth))
    )
  }

  # Vertices
  est_vert <- anomaly$mean
  # Observation locations
  true_obs <- as.matrix(Phi %*% true_vert)
  est_obs <- as.matrix(Phi %*% est_vert)
  naive_obs <- Y_mat - rowMeans(Y_mat)

  if (plot_results) {
    n_obs <- length(true_obs)
    # Thin to at most 20000 points per panel for fast plotting (deterministic, no RNG)
    keep <- unique(round(seq(1, n_obs, length.out = min(n_obs, 20000L))))
    scatter_df <- rbind(
      data.frame(method = "Model (posterior mean)", truth = true_obs[keep], estimate = est_obs[keep]),
      data.frame(method = "Naive baseline", truth = true_obs[keep], estimate = naive_obs[keep])
    )
    lim <- range(scatter_df$truth, scatter_df$estimate)

    p <- ggplot2::ggplot(scatter_df, ggplot2::aes(x = truth, y = estimate)) +
      ggplot2::geom_point(alpha = 0.15, size = 0.4) +
      ggplot2::geom_abline(slope = 1, intercept = 0, colour = "firebrick", linewidth = 0.8) +
      ggplot2::facet_wrap(~method, ncol = 2) +
      ggplot2::coord_equal(xlim = lim, ylim = lim) +
      ggplot2::labs(
        title = "Anomaly at the observation locations",
        subtitle = "Red: estimate = truth",
        x = "True anomaly",
        y = "Estimate"
      ) +
      diag_theme()
    print(p)
  }

  data.frame(
    Metric = c("RMSE", "MAE", "Correlation with truth"),
    Model_vertices = round(metrics(est_vert, true_vert), 4),
    Model_obs = round(metrics(est_obs, true_obs), 4),
    Naive_obs = round(metrics(naive_obs, true_obs), 4)
  )
}

#' Plot 2x2 Spatial Maps of the Anomaly for a Selected Year (ggplot2)
#'
#' Vertex values are projected onto a regular grid with the P1 finite element
#' basis (as in `plot_spde_field()`). The maps cover the bounding box of the
#' observation locations (`sim_data$locs`), or the whole mesh if the data are at
#' the vertices: vertices in the mesh extension have no data, so their large
#' posterior sd would otherwise dominate the colour scales.
#' The true anomaly and the posterior mean share one colour scale (limits = range
#' of both maps), so their colours can be compared directly. The posterior sd
#' (sequential palette) and the error (diverging palette centred at 0) each have
#' their own scale. Requires the ggplot2 and patchwork packages.
#'
#' @param mesh fmesher 2D mesh object
#' @param sim_data Output list of the toy data generator
#' @param anomaly Output of `summarize_anomaly(res)`
#' @param year Integer index of year to display
#' @param resolution Number of grid points per axis
#' @return Invisibly returns the patchwork object (it is also printed)
plot_spatial_year_diagnostics <- function(
    mesh,
    sim_data,
    anomaly,
    year = 1L,
    resolution = 150L
) {
  true_vert <- get_true_anomaly(sim_data)
  err <- anomaly$mean[, year] - true_vert[, year]

  # Regular grid over the data area, projected with the P1 basis
  area <- if (is.null(sim_data$locs)) mesh$loc else as.matrix(sim_data$locs)
  grid <- expand.grid(
    x = seq(min(area[, 1]), max(area[, 1]), length.out = resolution),
    y = seq(min(area[, 2]), max(area[, 2]), length.out = resolution)
  )
  A <- fmesher::fm_basis(mesh, loc = as.matrix(grid))
  inside <- Matrix::rowSums(A) > 0
  grid <- grid[inside, ]
  A <- A[inside, , drop = FALSE]

  on_grid <- function(values) as.vector(A %*% values)

  make_panel <- function(values, title, fill_scale) {
    df <- grid
    df$value <- on_grid(values)
    ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, fill = value)) +
      ggplot2::geom_raster() +
      ggplot2::coord_equal(expand = FALSE) +
      fill_scale +
      ggplot2::labs(title = title, x = NULL, y = NULL) +
      diag_theme(base_size = 10) +
      ggplot2::theme(
        legend.position = "right",
        legend.key.width = ggplot2::unit(0.3, "cm"),
        legend.key.height = ggplot2::unit(0.55, "cm"),
        legend.title = ggplot2::element_text(size = 9),
        legend.text = ggplot2::element_text(size = 8),
        plot.margin = ggplot2::margin(2, 4, 2, 4),
        panel.grid = ggplot2::element_blank()
      )
  }

  # Truth and posterior mean share one colour scale
  shared_lim <- range(on_grid(true_vert[, year]), on_grid(anomaly$mean[, year]))
  anomaly_scale <- function() {
    ggplot2::scale_fill_distiller(
      palette = "Spectral",
      direction = -1,
      limits = shared_lim,
      name = "Anomaly"
    )
  }
  err_max <- max(abs(on_grid(err)))

  p_true <- make_panel(true_vert[, year], "True anomaly", anomaly_scale())
  p_mean <- make_panel(anomaly$mean[, year], "Posterior mean", anomaly_scale())
  p_sd <- make_panel(
    anomaly$sd[, year],
    "Posterior sd",
    ggplot2::scale_fill_viridis_c(option = "magma", direction = -1, name = "sd")
  )
  p_err <- make_panel(
    err,
    "Error (mean - truth)",
    ggplot2::scale_fill_distiller(
      palette = "RdBu",
      direction = -1,
      limits = c(-err_max, err_max),
      name = "Error"
    )
  )

  # Top row: one shared legend for both maps; bottom row: one legend per map
  top_row <- patchwork::wrap_plots(p_true, p_mean, ncol = 2, guides = "collect")
  bottom_row <- patchwork::wrap_plots(p_sd, p_err, ncol = 2)
  p <- patchwork::wrap_plots(top_row, bottom_row, ncol = 1) +
    patchwork::plot_annotation(
      title = sprintf("Anomaly field, year %d", year),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 13))
    )
  print(p)
  invisible(p)
}

#' Run Complete Diagnostic Suite and Print Summary Tables
#'
#' @param res Output list returned by `run_gibbs_sampler()`
#' @param mesh fmesher 2D mesh object
#' @param sim_data Output list of the toy data generator
#' @param truth List of true values (kappa, sigma_beta, v1, v2, sigma_e)
#' @param num_burn_in Number of burn-in iterations
#' @param elapsed Optional elapsed time in seconds
#' @param year Year index for the spatial maps
#' @param vertex Vertex index for the beta_tilde and mu traces
#' @return Invisibly returns list of summary tables
run_all_diagnostics <- function(
    res,
    mesh,
    sim_data,
    truth,
    num_burn_in,
    elapsed = NULL,
    year = 1L,
    vertex = 1L
) {
  cat("\n=== MCMC RUN SUMMARY ===\n")
  if (!is.null(elapsed)) {
    cat(sprintf("Elapsed Time: %.2f s | ", elapsed))
  }
  cat(sprintf("Acceptance Rate (post burn-in): %.2f%%\n\n", res$acceptance_rate * 100))

  plot_mcmc_traces(res, sim_data, truth, num_burn_in, vertex = vertex)

  ess_tab <- plot_mcmc_acf_ess(res, num_burn_in, vertex = vertex)
  cat("=== EFFECTIVE SAMPLE SIZE (POST BURN-IN) ===\n")
  print(ess_tab, row.names = FALSE)
  cat("\n")

  cor_mat <- compute_hyper_correlations(res, num_burn_in)
  cat("=== POSTERIOR CORRELATIONS ===\n")
  print(cor_mat)
  cat("\n")

  est_tab <- summarize_hyperparameters(res, sim_data, truth, num_burn_in, vertex = vertex)
  cat("=== POSTERIOR ESTIMATES VS GROUND TRUTH ===\n")
  print(est_tab)
  cat("\n")

  anomaly <- summarize_anomaly(res)
  field_tab <- evaluate_field_accuracy(res, sim_data, anomaly, plot_results = TRUE)
  cat("=== ANOMALY FIELD RECONSTRUCTION ACCURACY ===\n")
  print(field_tab, row.names = FALSE)

  mu_rmse <- sqrt(mean((rowMeans(res$mu_samples) - sim_data$mu_true)^2))
  cat(sprintf("\nRMSE of posterior mean of mu (vertices): %.4f\n\n", mu_rmse))

  plot_spatial_year_diagnostics(mesh, sim_data, anomaly, year = year)

  invisible(list(
    ess = ess_tab,
    correlations = cor_mat,
    estimates = est_tab,
    field_metrics = field_tab,
    mu_rmse = mu_rmse,
    anomaly = anomaly
  ))
}
