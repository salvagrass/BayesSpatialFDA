#' Generate Synthetic Spatial Data with Year-to-Year Persistence (AR(1) anomalies)
#'
#' Same data-generating process as [generate_toy_data()], except that the yearly
#' anomalies are not independent: they follow a first-order autoregressive
#' process in time,
#' \deqn{\tilde{\beta}_1 = \xi_1, \qquad
#'       \tilde{\beta}_t = \rho \, \tilde{\beta}_{t-1} + \sqrt{1 - \rho^2} \, \xi_t,
#'       \quad t = 2, \dots, M,}
#' where \eqn{\xi_1, \dots, \xi_M} are independent SPDE fields
#' \eqn{\mathcal{N}(0, Q(\theta)^{-1})}, drawn exactly as in [generate_toy_data()].
#' The factor \eqn{\sqrt{1 - \rho^2}} keeps every year's spatial distribution equal
#' to \eqn{\mathcal{N}(0, Q(\theta)^{-1})}, so only the dependence between years
#' changes: \eqn{\mathrm{Cov}(\tilde{\beta}_t, \tilde{\beta}_{t+k}) = \rho^{k} Q(\theta)^{-1}}.
#' With `rho = 0` the output equals that of [generate_toy_data()] with the same
#' seed.
#'
#' @param mesh An `fm_mesh_2d` object representing the spatial domain.
#' @param num_year Integer. The number of years (`M`) to simulate. Default is `10`.
#' @param mu_prior_mean Numeric scalar or vector. Mean of the prior for `mu`. Default is `0.0`.
#' @param Q_mu Sparse matrix or scalar. The precision matrix for `mu`. Default is `1.0`.
#' @param kappa Numeric scalar. SPDE spatial scale parameter. Default is `1.0`.
#' @param sigma_beta Numeric scalar. Marginal standard deviation of the latent field. Default is `1.0`.
#' @param v1 Numeric scalar. First component of the anisotropy vector `v`. Default is `0.0`.
#' @param v2 Numeric scalar. Second component of the anisotropy vector `v`. Default is `0.0`.
#' @param rho Numeric scalar in `(-1, 1)`. Correlation between the anomalies of
#'   consecutive years at the same location. `0` gives independent years (the
#'   model's assumption); positive values give persistence. Default is `0.5`.
#' @param sigma_e Numeric scalar. Standard deviation of observation noise. Default is `0.1`.
#' @param locs Optional numeric matrix (`N_obs x 2` or `N_obs x 3`) of observation
#'   coordinates. If `NULL` (default), observations are generated directly at the mesh vertices.
#' @param seed Optional integer for reproducibility.
#'
#' @return A list containing `y_mat` (`N_obs x num_year`), `mu_true` (`V`),
#'   `beta_mat` (`V x num_year`), `Phi` (`N_obs x V`), `locs`, and `params`.
#' @export
generate_toy_data_ar1 <- function(
    mesh,
    num_year = 10,
    mu_prior_mean = 0.0,
    Q_mu = 1.0,
    kappa = 1.0,
    sigma_beta = 1.0,
    v1 = 0.0,
    v2 = 0.0,
    rho = 0.5,
    sigma_e = 0.1,
    locs = NULL,
    seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (length(rho) != 1L || !is.finite(rho) || abs(rho) >= 1) {
    stop("Arguments error: 'rho' must be a single number in (-1, 1).")
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

  # 3. Latent fields at mesh vertices (V x num_year), with AR(1) anomalies
  # Independent SPDE innovations, drawn as in generate_toy_data()
  xi_mat <- replicate(
    num_year,
    sample_normal_vector(
      mean = 0,
      mat = Q_beta,
      N = V,
      is_precision = TRUE
    )
  )
  # AR(1) recursion in time; every year keeps the marginal N(0, Q_beta^-1)
  anomaly_mat <- xi_mat
  if (num_year > 1) {
    for (t in 2:num_year) {
      anomaly_mat[, t] <- rho * anomaly_mat[, t - 1] + sqrt(1 - rho^2) * xi_mat[, t]
    }
  }
  beta_mat <- mu_true + anomaly_mat

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
      rho = rho,
      sigma_e = sigma_e,
      seed = seed
    )
  ))
}
