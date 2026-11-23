#' Generate Synthetic Spatial Data with Spatially Varying kappa and sigma_beta (non-stationarity)
#'
#' Same data-generating process as [generate_toy_data()], except that the SPDE
#' parameters `kappa` (inverse range) and `sigma_beta` (marginal standard
#' deviation) may vary over the domain. The anisotropy `v` stays constant.
#'
#' The yearly fields are drawn from the non-stationary SPDE precision
#' \deqn{Q = D \, (\tilde{C} K^2 + G_{H_v}) \, \tilde{C}^{-1} (\tilde{C} K^2 + G_{H_v}) \, D,}
#' with \eqn{K = \mathrm{diag}(\kappa_j)} and
#' \eqn{D = \mathrm{diag}\big(1 / (\sqrt{4\pi}\,\kappa_j \sigma_{\beta,j})\big)},
#' where \eqn{\tilde{C}} is the lumped diagonal mass matrix and \eqn{G_{H_v}} the
#' anisotropic stiffness matrix, assembled exactly as in [build_precision_matrix()].
#' With constant `kappa` and `sigma_beta` this is the precision of
#' [generate_toy_data()], and the output equals that of [generate_toy_data()]
#' with the same seed, so the spatial variation of the parameters is the only
#' difference between the two generators. When `kappa` varies slowly compared to
#' the range, the local range is `sqrt(8) / kappa(s)` and the local marginal
#' standard deviation is `sigma_beta(s)` (away from the mesh boundary).
#'
#' @param mesh An `fm_mesh_2d` object representing the spatial domain.
#' @param num_year Integer. The number of years (`M`) to simulate. Default is `10`.
#' @param mu_prior_mean Numeric scalar or vector. Mean of the prior for `mu`. Default is `0.0`.
#' @param Q_mu Sparse matrix or scalar. The precision matrix for `mu`. Default is `1.0`.
#' @param kappa SPDE spatial scale parameter: a positive scalar (stationary), a
#'   vector of length `V` (one value per mesh vertex), or a `function(loc)` that
#'   receives the vertex coordinates `mesh$loc` (`V x 3`) and returns `V` values.
#'   Default is `1.0`.
#' @param sigma_beta Marginal standard deviation of the latent field: scalar,
#'   vector of length `V`, or `function(loc)`, as for `kappa`. Default is `1.0`.
#' @param v1 Numeric scalar. First component of the anisotropy vector `v`. Default is `0.0`.
#' @param v2 Numeric scalar. Second component of the anisotropy vector `v`. Default is `0.0`.
#' @param sigma_e Numeric scalar. Standard deviation of observation noise. Default is `0.1`.
#' @param locs Optional numeric matrix (`N_obs x 2` or `N_obs x 3`) of observation
#'   coordinates. If `NULL` (default), observations are generated directly at the mesh vertices.
#' @param seed Optional integer for reproducibility.
#'
#' @return A list containing `y_mat` (`N_obs x num_year`), `mu_true` (`V`),
#'   `beta_mat` (`V x num_year`), `Phi` (`N_obs x V`), `locs`, and `params`.
#'   `params$kappa` and `params$sigma_beta` hold the values at the mesh vertices.
#' @export
generate_toy_data_nonstationary <- function(
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

  # Parameter values at the mesh vertices: a scalar, a vector of length V,
  # or a function of the vertex coordinates, strictly positive everywhere
  at_vertices <- function(x, name) {
    values <- as.numeric(if (is.function(x)) x(mesh$loc) else x)
    if (length(values) == 1L) {
      values <- rep(values, V)
    }
    if (length(values) != V) {
      stop(sprintf(
        "Arguments error: '%s' must have length 1 or %d (number of vertices), got %d.",
        name, V, length(values)
      ))
    }
    if (anyNA(values) || any(values <= 0)) {
      stop(sprintf("Arguments error: '%s' must be strictly positive at every vertex.", name))
    }
    values
  }
  kappa_vec <- at_vertices(kappa, "kappa")
  sigma_vec <- at_vertices(sigma_beta, "sigma_beta")

  # 1. Climatological mean field on mesh vertices
  mu_true <- sample_normal_vector(
    mean = mu_prior_mean,
    mat = Q_mu,
    N = V,
    is_precision = TRUE
  )

  # 2. Non-stationary SPDE precision matrix for yearly fields
  # (same assembly as build_precision_matrix(), with kappa and sigma_beta per vertex)
  matrices <- build_spde_matrices(mesh, manifold = "planar")
  H <- anisotropy_from_v(v1, v2)

  global_K <- H$H11 *
    matrices$K11 +
    H$H12 * matrices$K12 +
    H$H22 * matrices$K22

  c_lumped <- rowSums(matrices$C)
  C_tilde_inv <- Matrix::Diagonal(x = 1 / c_lumped)

  core_mat <- global_K + Matrix::Diagonal(x = kappa_vec^2 * c_lumped)
  D <- Matrix::Diagonal(x = 1 / (sqrt(4 * pi) * kappa_vec * sigma_vec))

  Q_beta <- D %*% core_mat %*% C_tilde_inv %*% core_mat %*% D

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
      kappa = kappa_vec,
      sigma_beta = sigma_vec,
      v1 = v1,
      v2 = v2,
      sigma_e = sigma_e,
      seed = seed
    )
  ))
}
