# experiments/gibbs_sampler/diagnostics.R
# Diagnostic plots and accuracy metrics for the Anisotropic SPDE Gibbs Sampler

#' Plot Full MCMC Traceplots with Burn-in and Ground Truth
#'
#' @param res Output list returned by `run_gibbs_sampler()`
#' @param truth List of true values (`sigma2`, `tau2`, `kappa2`, `gamma`, `theta`, `beta1`)
plot_mcmc_traces <- function(res, truth = list()) {
  burnin <- res$diagnostics$burnin
  d <- res$draws

  trace_list <- list(
    list(x = d$sigma2, main = "sigma^2 (Noise)", true = truth$sigma2),
    list(x = d$tau2, main = "tau^2", true = truth$tau2),
    list(x = d$kappa2, main = "kappa^2", true = truth$kappa2),
    list(x = d$gamma, main = "gamma (Anisotropy Ratio)", true = truth$gamma),
    list(x = d$theta, main = "theta (Angle)", true = truth$theta),
    list(x = d$beta1, main = "beta[1,1] (Anomaly)", true = truth$beta1)
  )

  old_par <- par(mfrow = c(2, 3))
  on.exit(par(old_par), add = TRUE)

  for (tr in trace_list) {
    plot(tr$x, type = "l", main = tr$main, xlab = "Iteration", ylab = "")
    abline(v = burnin, lty = 2, col = "grey40")
    if (!is.null(tr$true)) {
      abline(h = tr$true, col = "red", lwd = 2)
    }
  }
}

#' Build Post Burn-In Matrix on Unconstrained/Log Scale
#'
#' @param res Output list returned by `run_gibbs_sampler()`
#' @return Numeric matrix of post burn-in draws
get_post_draws_matrix <- function(res) {
  idx <- res$diagnostics$idx_post
  d <- res$draws

  cbind(
    "log sigma^2" = log(d$sigma2[idx]),
    "log tau^2" = log(d$tau2[idx]),
    "log kappa^2" = log(d$kappa2[idx]),
    "log gamma" = log(d$gamma[idx]),
    "theta" = d$theta[idx],
    "beta[1,1]" = d$beta1[idx]
  )
}

#' Plot ACF and Compute Effective Sample Size (ESS)
#'
#' @param res Output list returned by `run_gibbs_sampler()`
#' @param lag_max Maximum lag for ACF plots
#' @return Data frame with ESS and total post-burn-in iterations
plot_mcmc_acf_ess <- function(res, lag_max = 60L) {
  post_draws <- get_post_draws_matrix(res)

  old_par <- par(mfrow = c(2, 3))
  on.exit(par(old_par), add = TRUE)

  for (j in seq_len(ncol(post_draws))) {
    acf(
      post_draws[, j],
      lag.max = lag_max,
      main = paste("ACF:", colnames(post_draws)[j])
    )
  }

  ess <- coda::effectiveSize(coda::mcmc(post_draws))
  ess_df <- data.frame(
    Parameter = names(ess),
    ESS = round(as.numeric(ess)),
    Out_of = nrow(post_draws),
    row.names = NULL
  )
  return(ess_df)
}

#' Plot Hyperparameter Pairs and Return Correlation Matrix
#'
#' @param res Output list returned by `run_gibbs_sampler()`
#' @return Rounded correlation matrix of the 5 hyperparameters
plot_hyper_correlations <- function(res) {
  post_draws <- get_post_draws_matrix(res)[, 1:5, drop = FALSE]
  pairs(post_draws, pch = ".", col = adjustcolor("black", 0.3))
  return(round(cor(post_draws), 2))
}

#' Summarize Posterior Estimates vs Ground Truth
#'
#' @param res Output list returned by `run_gibbs_sampler()`
#' @param truth List of true values (`sigma2`, `tau2`, `kappa2`, `gamma`, `theta`, `beta1`)
#' @return Matrix with True, Mean, 2.5%, Median, and 97.5% quantiles
summarize_hyperparameters <- function(res, truth = list()) {
  idx <- res$diagnostics$idx_post
  d <- res$draws

  summ_vec <- function(x, true_val) {
    q <- quantile(x[idx], c(0.025, 0.5, 0.975))
    c(
      True = if (is.null(true_val)) NA_real_ else true_val,
      Mean = mean(x[idx]),
      `2.5%` = q[[1]],
      Median = q[[2]],
      `97.5%` = q[[3]]
    )
  }

  est_tab <- rbind(
    "sigma^2" = summ_vec(d$sigma2, truth$sigma2),
    "tau^2" = summ_vec(d$tau2, truth$tau2),
    "kappa^2" = summ_vec(d$kappa2, truth$kappa2),
    "gamma" = summ_vec(d$gamma, truth$gamma),
    "theta" = summ_vec(d$theta, truth$theta),
    "beta_11" = summ_vec(d$beta1, truth$beta1)
  )
  return(round(est_tab, 4))
}

#' Evaluate Field Reconstruction Accuracy vs Naive Baseline
#'
#' Computes RMSE, MAE, Correlation, plots Model vs Naive scatter and RMSE by year.
#'
#' @param res Output list returned by `run_gibbs_sampler()`
#' @param Y_mat Observation matrix (N x M)
#' @param true_anomaly_mat True zero-mean anomaly matrix (N x M)
#' @param plot_results Logical; whether to generate scatter and yearly RMSE barplot
#' @return Data frame comparing Model vs Naive baseline metrics
evaluate_field_accuracy <- function(
  res,
  Y_mat,
  true_anomaly_mat,
  plot_results = TRUE
) {
  beta_mean <- res$summary$beta_mean
  err <- beta_mean - true_anomaly_mat

  rmse <- sqrt(mean(err^2))
  mae <- mean(abs(err))
  corr <- cor(as.vector(beta_mean), as.vector(true_anomaly_mat))

  # Naive baseline: each year minus node's temporal average
  naive <- Y_mat - rowMeans(Y_mat)
  rmse_naive <- sqrt(mean((naive - true_anomaly_mat)^2))
  mae_naive <- mean(abs(naive - true_anomaly_mat))
  corr_naive <- cor(as.vector(naive), as.vector(true_anomaly_mat))

  if (plot_results) {
    old_par <- par(mfrow = c(1, 2))
    plot(
      true_anomaly_mat,
      beta_mean,
      pch = 16,
      cex = 0.3,
      col = adjustcolor("black", 0.3),
      xlab = "True beta",
      ylab = "Posterior mean",
      main = "Model vs truth"
    )
    abline(0, 1, col = "red", lwd = 2)

    plot(
      true_anomaly_mat,
      naive,
      pch = 16,
      cex = 0.3,
      col = adjustcolor("black", 0.3),
      xlab = "True beta",
      ylab = "Naive estimate",
      main = "Naive baseline vs truth"
    )
    abline(0, 1, col = "red", lwd = 2)
    par(old_par)

    M <- ncol(Y_mat)
    rmse_year <- sqrt(colMeans(err^2))
    barplot(
      rmse_year,
      names.arg = seq_len(M),
      xlab = "Year",
      ylab = "RMSE",
      main = "RMSE of posterior mean by year",
      col = "grey70",
      border = NA
    )
    abline(h = rmse, col = "red", lwd = 2, lty = 2)
  }

  metrics_df <- data.frame(
    Metric = c("RMSE", "MAE", "Correlation with truth"),
    Model = round(c(rmse, mae, corr), 4),
    Naive_baseline = round(c(rmse_naive, mae_naive, corr_naive), 4)
  )
  return(metrics_df)
}

#' Plot 2x2 Spatial Maps for a Selected Year
#'
#' @param mesh fmesher 2D mesh object
#' @param res Output list returned by `run_gibbs_sampler()`
#' @param true_anomaly_mat True zero-mean anomaly matrix (N x M)
#' @param year Integer index of year to display
plot_spatial_year_diagnostics <- function(
  mesh,
  res,
  true_anomaly_mat,
  year = 1L
) {
  beta_mean <- res$summary$beta_mean
  beta_sd <- res$summary$beta_sd
  err <- beta_mean - true_anomaly_mat

  old_par <- par(mfrow = c(2, 2))
  on.exit(par(old_par), add = TRUE)

  plot_spde_field(
    mesh,
    true_anomaly_mat[, year],
    main_title = paste("True anomaly, year", year)
  )
  plot_spde_field(mesh, beta_mean[, year], main_title = "Posterior mean")
  plot_spde_field(mesh, beta_sd[, year], main_title = "Posterior sd")
  plot_spde_field(
    mesh,
    err[, year],
    main_title = "Error (posterior mean - truth)"
  )
}

#' Run Complete Diagnostic Suite and Print Summary Tables
#'
#' @param res Output list returned by `run_gibbs_sampler()`
#' @param mesh fmesher 2D mesh object
#' @param Y_mat Observation matrix (N x M)
#' @param true_anomaly_mat True zero-mean anomaly matrix (N x M)
#' @param truth List of true hyperparameter values
#' @param year Year index for spatial map panel
#' @return Invisibly returns list of summary tables (ess, cor, estimates, field_metrics)
run_all_diagnostics <- function(
  res,
  mesh,
  Y_mat,
  true_anomaly_mat,
  truth = list(),
  year = 1L
) {
  cat(sprintf(
    "\n=== MCMC RUN SUMMARY ===\nElapsed Time: %.2f s | Acceptance Rate: %.2f%%\n\n",
    res$diagnostics$elapsed_sec,
    res$diagnostics$acc_rate * 100
  ))

  plot_mcmc_traces(res, truth = truth)

  ess_tab <- plot_mcmc_acf_ess(res)
  cat("=== EFFECTIVE SAMPLE SIZE (POST BURN-IN) ===\n")
  print(ess_tab, row.names = FALSE)
  cat("\n")

  cor_mat <- plot_hyper_correlations(res)
  cat("=== POSTERIOR CORRELATIONS ===\n")
  print(cor_mat)
  cat("\n")

  est_tab <- summarize_hyperparameters(res, truth = truth)
  cat("=== POSTERIOR ESTIMATES VS GROUND TRUTH ===\n")
  print(est_tab)
  cat("\n")

  field_tab <- evaluate_field_accuracy(
    res,
    Y_mat,
    true_anomaly_mat,
    plot_results = TRUE
  )
  cat("=== FIELD RECONSTRUCTION ACCURACY ===\n")
  print(field_tab, row.names = FALSE)
  cat("\n")

  plot_spatial_year_diagnostics(mesh, res, true_anomaly_mat, year = year)

  invisible(list(
    ess = ess_tab,
    correlations = cor_mat,
    estimates = est_tab,
    field_metrics = field_tab
  ))
}
