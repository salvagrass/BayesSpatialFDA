#' Generate Synthetic Spatial Data with Matern Smoothness nu (misspecified smoothness)
#'
#' Same data-generating process as [generate_toy_data()], except that the yearly
#' latent fields are drawn from a dense anisotropic Matern covariance with
#' smoothness `nu` instead of the SPDE precision (which fixes `nu = 1`).
#'
#' The correlation range `sqrt(8) / kappa` and the anisotropy `H_v` keep the same
#' meaning as in the SPDE model:
#' * the Matern correlation is evaluated in the metric `d(h) = sqrt(h' H_v^-1 h)`,
#'   which is the geometry induced by the operator `kappa^2 - div(H_v grad)`;
#' * the inverse range is rescaled to `kappa_nu = sqrt(nu) * kappa`, so that
#'   `sqrt(8 * nu) / kappa_nu = sqrt(8) / kappa` for every `nu`.
#'
#' `nu = 1` gives the continuous-domain counterpart of the SPDE field (without
#' finite element discretization and boundary effects). The covariance is dense,
#' so keep the mesh to a few thousand vertices.
#'
#' @param mesh An `fm_mesh_2d` object representing the spatial domain.
#' @param num_year Integer. The number of years (`M`) to simulate. Default is `10`.
#' @param mu_prior_mean Numeric scalar or vector. Mean of the prior for `mu`. Default is `0.0`.
#' @param Q_mu Sparse matrix or scalar. The precision matrix for `mu`. Default is `1.0`.
#' @param kappa Numeric scalar. Spatial scale parameter; the correlation range is
#'   `sqrt(8) / kappa` for every `nu`. Default is `1.0`.
#' @param sigma_beta Numeric scalar. Marginal standard deviation of the latent field. Default is `1.0`.
#' @param v1 Numeric scalar. First component of the anisotropy vector `v`. Default is `0.0`.
#' @param v2 Numeric scalar. Second component of the anisotropy vector `v`. Default is `0.0`.
#' @param nu Numeric scalar. Matern smoothness. Default is `0.5` (exponential covariance).
#' @param sigma_e Numeric scalar. Standard deviation of observation noise. Default is `0.1`.
#' @param locs Optional numeric matrix (`N_obs x 2` or `N_obs x 3`) of observation
#'   coordinates. If `NULL` (default), observations are generated directly at the mesh vertices.
#' @param seed Optional integer for reproducibility.
#'
#' @return A list containing `y_mat` (`N_obs x num_year`), `mu_true` (`V`),
#'   `beta_mat` (`V x num_year`), `Phi` (`N_obs x V`), `locs`, and `params`.
#' @export
generate_toy_data_matern <- function(
    mesh,
    num_year = 10,
    mu_prior_mean = 0.0,
    Q_mu = 1.0,
    kappa = 1.0,
    sigma_beta = 1.0,
    v1 = 0.0,
    v2 = 0.0,
    nu = 0.5,
    sigma_e = 0.1,
    locs = NULL,
    seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (nu <= 0) {
    stop("Arguments error: the Matern smoothness 'nu' must be positive.")
  }

  V <- nrow(mesh$loc)

  # 1. Climatological mean field on mesh vertices
  mu_true <- sample_normal_vector(
    mean = mu_prior_mean,
    mat = Q_mu,
    N = V,
    is_precision = TRUE
  )

  # 2. Dense anisotropic Matern covariance for yearly fields
  # Since det(H_v) = 1, H_v^-1 = [H22, -H12; -H12, H11], and the squared
  # distance in the metric of H_v^-1 is h' H_v^-1 h.
  H <- anisotropy_from_v(v1, v2)
  dx <- outer(mesh$loc[, 1], mesh$loc[, 1], "-")
  dy <- outer(mesh$loc[, 2], mesh$loc[, 2], "-")
  dist_mat <- sqrt(pmax(H$H22 * dx^2 - 2 * H$H12 * dx * dy + H$H11 * dy^2, 0))

  scaled_d <- sqrt(nu) * kappa * dist_mat
  Sigma_beta <- sigma_beta^2 *
    (2^(1 - nu) / gamma(nu)) *
    scaled_d^nu *
    besselK(scaled_d, nu)
  diag(Sigma_beta) <- sigma_beta^2
  # Small nugget: very smooth fields (large nu) are numerically singular
  diag(Sigma_beta) <- diag(Sigma_beta) + 1e-8 * sigma_beta^2

  # 3. Latent fields at mesh vertices (V x num_year)
  beta_mat <- replicate(
    num_year,
    sample_normal_vector(
      mean = mu_true,
      mat = Sigma_beta,
      N = V,
      is_precision = FALSE
    )
  )

  # 4. Projection to observation locations (if provided)
  if (is.null(locs)) {
    Phi <- Matrix::Diagonal(V)
    signal_mat <- beta_mat
  } else {
    Phi <- fmesher::fm_basis(mesh, loc = locs)
    stopifnot(abs(rowSums(Phi) - 1) < 1e-12)
    signal_mat <- as.matrix(Phi %*% beta_mat)
  }

  # 5. Observation noise and final Y (N_obs x num_year)
  N_obs <- nrow(signal_mat)
  epsilon <- matrix(
    data = rnorm(N_obs * num_year, 0, sigma_e),
    nrow = N_obs,
    ncol = num_year
  )
  y_mat <- signal_mat + epsilon

  return(list(
    y_mat = y_mat,
    mu_true = mu_true,
    beta_mat = beta_mat,
    Phi = Phi,
    locs = locs,
    params = list(
      kappa = kappa,
      sigma_beta = sigma_beta,
      v1 = v1,
      v2 = v2,
      nu = nu,
      sigma_e = sigma_e,
      seed = seed
    )
  ))
}
