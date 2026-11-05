#' @title Generate Synthetic Spatial Data for Multi-Year SPDE Models
#' @description
#' This function generates synthetic spatial data on a 2D mesh domain for multiple years,
#' following the exact model formulation:
#' Y_i = a_i + mu + beta_tilde_i + epsilon_i
#'
#' Where for each year i = 1, ..., num_year:
#' - a_i is the global temporal effect for year i, drawn from Normal(a_prior_mean, a_prior_sd^2).
#' - mu is the long-term average temperature field (constant across years), drawn from Normal(mu_prior_mean, Q_mu^-1).
#' - beta_tilde_i is the zero-mean spatially correlated anomaly for year i, drawn via SPDE precision Q(tau, h1, h2).
#' - epsilon_i is the unstructured measurement error for year i.
#'
#' @param mesh An `fm_mesh_2d` object representing the spatial domain.
#' @param num_year Integer. The number of years (M) to simulate. Default is 10.
#' @param a_prior_mean Numeric scalar. Mean of the prior for the global temporal effects a_i. Default is 0.0.
#' @param a_prior_sd Numeric scalar. Standard deviation of the prior for a_i. Default is 10.0.
#' @param mu_prior_mean Numeric scalar or vector. Mean of the independent prior for the mu field. Default is 0.0.
#' @param Q_mu Sparse matrix. The precision matrix for mu. Default is 1.0 (Scalar values imply diagonal matrices).
#' @param tau Numeric scalar. SPDE precision parameter. Default is 1.0.
#' @param k Numeric scalar. SPDE geometric reaction parameter. Default is 1.0.
#' @param sigma_e Numeric scalar. Standard deviation of observation noise. Default is 0.1.
#' @param seed Optional integer for reproducibility.
#'
#' @return A list containing the observations and the exact ground truth components:
#' \itemize{
#'   \item \code{y_mat}: A matrix of observations (num_year x N).
#'   \item \code{a_true}: A vector of the global temporal effects (length num_year).
#'   \item \code{mu_true}: A vector of the long-term spatial mean field (length N).
#'   \item \code{beta_tilde_mat}: A matrix of the zero-mean spatial anomalies (num_year x N).
#'   \item \code{latent_field_mat}: A matrix of the true latent field without noise (num_year x N).
#' }
#' @export
generate_toy_data <- function(
  mesh,
  num_year = 10,
  a_prior_mean = 0.0,
  a_prior_sd = 10.0,
  mu_prior_mean = 0.0,
  Q_mu = 1.0,
  tau = 1.0,
  k = 1.0,
  sigma_e = 0.1,
  seed = NULL
) {
  # Setting seed if given
  if (!is.null(seed)) {
    set.seed(seed)
  }

  N <- nrow(mesh$loc)

  # 1. Global temporal effect vector (vector)
  a_true <- sample_normal_vector(
    mean = a_prior_mean,
    mat = a_prior_sd,
    N = num_year
  )

  mu_true <- sample_normal_vector(
    mean = mu_prior_mean,
    mat = Q_mu,
    N = N,
    is_precision = TRUE
  )

  # 3. SPDE Precision for beta_tilde_i
  matrices <- build_spde_matrices(mesh, manifold = "planar")

  core_mat <- matrices$K + matrices$C
  diag_C <- diag(matrices$C)
  Q_beta <- (tau^2) * core_mat %*% Diagonal(x = 1 / diag_C) %*% core_mat

  beta_tilde_mat <- replicate(
    num_year,
    sample_normal_vector(
      mean = 0.0,
      mat = Q_beta,
      N = N,
      is_precision = TRUE
    )
  )

  # 4. Total Latent field
  latent_field_mat <- outer(mu_true, a_true, "+") + beta_tilde_mat

  # 5. Observation noise and final Y
  epsilon <- matrix(
    data = rnorm(N * num_year, 0, sigma_e),
    nrow = N,
    ncol = num_year
  )
  y_mat <- latent_field_mat + epsilon

  return(list(
    y_mat = y_mat,
    a_true = a_true,
    mu_true = mu_true,
    beta_tilde_mat = beta_tilde_mat,
    latent_field_mat = latent_field_mat
  ))
}
