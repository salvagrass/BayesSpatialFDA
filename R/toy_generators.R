#' Generate Synthetic Spatial Data for Multi-Year SPDE Models
#'
#' @param mesh An `fm_mesh_2d` object representing the spatial domain.
#' @param num_year Integer. The number of years (`M`) to simulate. Default is `10`.
#' @param mu_prior_mean Numeric scalar or vector. Mean of the prior for `mu`. Default is `0.0`.
#' @param Q_mu Sparse matrix or scalar. The precision matrix for `mu`. Default is `1.0`.
#' @param kappa Numeric scalar. SPDE spatial scale parameter. Default is `1.0`.
#' @param sigma_beta Numeric scalar. Marginal standard deviation of the latent field. Default is `1.0`.
#' @param v1 Numeric scalar. First component of the anisotropy vector `v`. Default is `0.0`.
#' @param v2 Numeric scalar. Second component of the anisotropy vector `v`. Default is `0.0`.
#' @param sigma_e Numeric scalar. Standard deviation of observation noise. Default is `0.1`.
#' @param locs Optional numeric matrix (`N_obs x 2` or `N_obs x 3`) of observation
#'   coordinates. If `NULL` (default), observations are generated directly at the mesh vertices.
#' @param seed Optional integer for reproducibility.
#'
#' @return A list containing `y_mat` (`N_obs x num_year`), `mu_true` (`V`),
#'   `beta_mat` (`V x num_year`), `Phi` (`N_obs x V`), `locs`, and `params`.
#' @export
generate_toy_data <- function(
  mesh,
  num_year = 10,
  mu_prior_mean = 0.0,
  Q_mu = 1.0,
  kappa = 1.0,
  sigma_beta = 1.0,
  v1 = 0.0,
  v2 = 0.0,
  sigma_e = 0.1,
  locs = NULL,
  seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  V <- nrow(mesh$loc)

  # 1. Climatological mean field on mesh vertices
  mu_true <- sample_normal_vector(
    mean = mu_prior_mean,
    mat = Q_mu,
    N = V,
    is_precision = TRUE
  )

  # 2. SPDE precision matrix for yearly fields
  matrices <- build_spde_matrices(mesh, manifold = "planar")
  Q_beta <- build_precision_matrix(matrices, kappa, sigma_beta, v1, v2)

  # 3. Latent fields at mesh vertices (V x num_year)
  beta_mat <- replicate(
    num_year,
    sample_normal_vector(
      mean = mu_true,
      mat = Q_beta,
      N = V,
      is_precision = TRUE
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
      sigma_e = sigma_e,
      seed = seed
    )
  ))
}
