#' @title Generate Synthetic Spatial Data for Multi-Year SPDE Models
#' @description
#' This function generates synthetic spatial data on a 2D mesh domain for multiple years,
#' following the identifiable model formulation:
#' Y_{t} = beta_{t} + epsilon_{t}
#'
#' Where for each year t = 1, ..., num_year:
#' - mu is the long-term average temperature field, drawn from Normal(mu_prior_mean, Q_mu^-1).
#' - beta_{t} is the spatially correlated latent field for year t, drawn from Normal(mu, Q_beta^-1).
#' - Q_beta is the anisotropic SPDE precision matrix governed by tau, k (kappa), gamma, and theta.
#' - epsilon_{t} is the unstructured measurement error (Gaussian white noise).
#'
#' @param mesh An `fm_mesh_2d` object representing the spatial domain.
#' @param num_year Integer. The number of years (M) to simulate. Default is 10.
#' @param mu_prior_mean Numeric scalar or vector. Mean of the independent prior for the mu field. Default is 0.0.
#' @param Q_mu Sparse matrix or scalar. The precision matrix for mu. Default is 1.0.
#' @param tau Numeric scalar. SPDE precision parameter. Default is 1.0.
#' @param k Numeric scalar. SPDE spatial scale parameter (kappa). Default is 1.0.
#' @param gamma Numeric scalar. Anisotropy ratio parameter. Default is 1.0 (Isotropic).
#' @param theta Numeric scalar. Anisotropy rotation angle in radians. Default is 0.0.
#' @param sigma_e Numeric scalar. Standard deviation of observation noise. Default is 0.1.
#' @param seed Optional integer for reproducibility.
#'
#' @return A list containing the observations and the exact ground truth components:
#' \itemize{
#'   \item \code{y_mat}: A matrix of observations (N x num_year).
#'   \item \code{mu_true}: A vector of the long-term spatial mean field (length N).
#'   \item \code{beta_mat}: A matrix of the true spatial latent fields (N x num_year).
#' }
#' @export
generate_toy_data <- function(
  mesh,
  num_year = 10,
  mu_prior_mean = 0.0,
  Q_mu = 1.0,
  tau = 1.0,
  k = 1.0,
  gamma = 1.0,
  theta = 0.0,
  sigma_e = 0.1,
  seed = NULL
) {
  # Setting seed if given
  if (!is.null(seed)) {
    set.seed(seed)
  }

  N <- nrow(mesh$loc)

  mu_true <- sample_normal_vector(
    mean = mu_prior_mean,
    mat = Q_mu,
    N = N,
    is_precision = TRUE
  )

  # 3. SPDE Precision for beta_tilde_i
  matrices <- build_spde_matrices(mesh, manifold = "planar")

  H11 <- gamma * (cos(theta)^2) + (1 / gamma) * (sin(theta)^2)
  H22 <- gamma * (sin(theta)^2) + (1 / gamma) * (cos(theta)^2)
  H12 <- (gamma - 1 / gamma) * sin(theta) * cos(theta)

  global_K <- H11 * matrices$K11 + H12 * matrices$K12 + H22 * matrices$K22
  core_mat <- global_K + (k^2) * matrices$C
  diag_C <- diag(matrices$C)
  Q_beta <- (tau^2) * core_mat %*% Diagonal(x = 1 / diag_C) %*% core_mat

  # 4. Total Latent field
  beta_mat <- replicate(
    num_year,
    sample_normal_vector(
      mean = mu_true,
      mat = Q_beta,
      N = N,
      is_precision = TRUE
    )
  )

  # 5. Observation noise and final Y
  epsilon <- matrix(
    data = rnorm(N * num_year, 0, sigma_e),
    nrow = N,
    ncol = num_year
  )
  y_mat <- beta_mat + epsilon

  return(list(
    y_mat = y_mat,
    mu_true = mu_true,
    beta_mat = beta_mat
  ))
}


# ==============================================================================
#' @title Generate Synthetic Spatial Data for Multi-Year SPDE Models (PC-prior parameterization)
#' @description
#' This function generates synthetic spatial data on a 2D mesh domain for multiple years,
#' following the identifiable model formulation:
#' Y_{t} = beta_{t} + epsilon_{t}
#'
#' Where for each year t = 1, ..., num_year:
#' - mu is the long-term average temperature field, drawn from Normal(mu_prior_mean, Q_mu^-1).
#' - beta_{t} is the spatially correlated latent field for year t, drawn from Normal(mu, Q_beta^-1).
#' - Q_beta is the anisotropic SPDE precision matrix governed by kappa, v = (v1, v2),
#'   and sigma_beta, following the half-angle parameterization of
#'   Llamazares-Elias, Latz & Lindgren:
#'   Q_beta = 1 / (4*pi*kappa^2*sigma_beta^2) * (kappa^2*C + G_H) C^-1 (kappa^2*C + G_H)
#' - epsilon_{t} is the unstructured measurement error (Gaussian white noise).
#'
#' @param mesh An `fm_mesh_2d` object representing the spatial domain.
#' @param num_year Integer. The number of years (M) to simulate. Default is 10.
#' @param mu_prior_mean Numeric scalar or vector. Mean of the independent prior for the mu field. Default is 0.0.
#' @param Q_mu Sparse matrix or scalar. The precision matrix for mu. Default is 1.0.
#' @param kappa Numeric scalar. SPDE spatial scale parameter; the correlation range is sqrt(8)/kappa. Default is 1.0.
#' @param v1 Numeric scalar. First component of the anisotropy vector v. Default is 0.0 (Isotropic).
#' @param v2 Numeric scalar. Second component of the anisotropy vector v. Default is 0.0 (Isotropic).
#' @param sigma_beta Numeric scalar. Marginal standard deviation of the latent field. Default is 1.0.
#' @param sigma_e Numeric scalar. Standard deviation of observation noise. Default is 0.1.
#' @param seed Optional integer for reproducibility.
#'
#' @return A list containing the observations and the exact ground truth components:
#' \itemize{
#'   \item \code{y_mat}: A matrix of observations (N x num_year).
#'   \item \code{mu_true}: A vector of the long-term spatial mean field (length N).
#'   \item \code{beta_mat}: A matrix of the true spatial latent fields (N x num_year).
#' }
#' @export
generate_toy_data_pc_prior <- function(
    mesh,
    num_year = 10,
    mu_prior_mean = 0.0,
    Q_mu = 1.0,
    kappa = 1.0,
    v1 = 0.0,
    v2 = 0.0,
    sigma_beta = 1.0,
    sigma_e = 0.1,
    seed = NULL
) {
  # Setting seed if given
  if (!is.null(seed)) {
    set.seed(seed)
  }

  N <- nrow(mesh$loc)

  mu_true <- sample_normal_vector(
    mean = mu_prior_mean,
    mat = Q_mu,
    N = N,
    is_precision = TRUE
  )

  # 3. SPDE Precision for beta_tilde_i
  matrices <- build_spde_matrices(mesh, manifold = "planar")

  # Anisotropy via the half-angle parameterization of Llamazares-Elias, Latz &
  # Lindgren: H_v = cosh(||v||) I + sinh(||v||)/||v|| * [v1 v2; v2 -v1].
  # This replaces the previous (gamma, theta) rotation/ratio parameterization.
  # The anisotropy strength is ||v|| (ratio a = exp(||v||))
  # and the direction is the half-angle arg(v)/2.
  norm_v <- sqrt(v1^2 + v2^2)
  if (norm_v < .Machine$double.eps) {
    # Isotropic limit: sinh(||v||)/||v|| -> 1 as ||v|| -> 0, H_v -> I
    sinhc_v <- 1.0
  } else {
    sinhc_v <- sinh(norm_v) / norm_v
  }
  h1 <- cosh(norm_v) + sinhc_v * v1
  h2 <- cosh(norm_v) - sinhc_v * v1
  h12 <- sinhc_v * v2

  global_K <- h1 * matrices$K11 + h12 * matrices$K12 + h2 * matrices$K22
  core_mat <- global_K + (kappa^2) * matrices$C
  diag_C <- diag(matrices$C)
  # Normalization 1/(4*pi*kappa^2*sigma_beta^2) fixes the marginal variance of
  # the field to sigma_beta^2.
  Q_beta <- (1 / (4 * pi * kappa^2 * sigma_beta^2)) *
    core_mat %*%
    Diagonal(x = 1 / diag_C) %*%
    core_mat

  # 4. Total Latent field
  beta_mat <- replicate(
    num_year,
    sample_normal_vector(
      mean = mu_true,
      mat = Q_beta,
      N = N,
      is_precision = TRUE
    )
  )

  # 5. Observation noise and final Y
  epsilon <- matrix(
    data = rnorm(N * num_year, 0, sigma_e),
    nrow = N,
    ncol = num_year
  )
  y_mat <- beta_mat + epsilon

  return(list(
    y_mat = y_mat,
    mu_true = mu_true,
    beta_mat = beta_mat
  ))
}

