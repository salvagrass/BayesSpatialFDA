#' Sample from a Multivariate Normal Distribution
#'
#' Generates a single random sample from a multivariate normal distribution.
#' It handles scalar inputs for both mean and covariance efficiently using sparse matrices.
#'
#' @param mean A numeric vector representing the mean. If a scalar is provided,
#'   it will be repeated `N` times.
#' @param mat A numeric matrix representing the covariance matrix or its inverse. If a scalar
#'   is provided, it is treated as the variance of independent variables and an
#'   internal sparse diagonal matrix is generated.
#' @param N An optional integer specifying the dimension of the vector. If `NULL`,
#'   it is inferred from the length of `mean` or the number of rows of `covariance`.
#'
#' @return A numeric vector of length `N` representing the generated sample.
#' @importFrom Matrix Diagonal Matrix
#' @export
sample_normal_vector <- function(mean, mat, N = NULL, is_precision = FALSE) {
  # Computing the dimensions of the vector if scalars are passed
  if (is.null(N)) {
    if (length(mean) > 1) {
      N <- length(mean)
    } else if (is.matrix(mat) || inherits(mat, "Matrix")) {
      N <- nrow(mat)
    } else {
      stop("Arguments error")
    }
  }

  if (length(mean) == 1) {
    mean_vec <- rep(mean, N)
  } else {
    mean_vec <- mean
  }

  if (length(mat) == 1) {
    mat_sym <- Diagonal(n = N, x = mat)
  } else {
    mat_sym <- Matrix(mat)
    if (!isSymmetric(mat_sym)) {
      stop("Arguments error: Covariance matrix must be symmetric.")
    }
  }

  # Add a tiny jitter on the diagonal for numerical stability
  mat_sym <- mat_sym + Matrix::Diagonal(N, 1e-9)

  # Sample generation, to generate samples with a non identiy covariance we use the Cholesky factorization
  Z <- rnorm(N)

  U <- tryCatch(
    {
      chol(mat_sym)
    },
    error = function(e) {
      stop("Error: ", e$message)
    }
  )

  # The mathematical branching
  if (is_precision) {
    # If mat is a Precision matrix Q: Q = U^T U.
    # Solving U * x = Z yields a vector with variance Q^(-1).
    samples <- as.vector(solve(U, Z)) + mean_vec
  } else {
    # If mat is a Covariance matrix Sigma: Sigma = U^T U.
    # Multiplying Z %*% U yields a vector with variance Sigma.
    samples <- as.vector(Z %*% U) + mean_vec
  }
  return(samples)
}


#' Plot a continuous SPDE spatial field using mesh and 1st-order FEM basis functions
#'
#' @description Interpolates and plots a continuous spatial map from mesh vertex values.
#' @param mesh The mesh object created with fmesher.
#' @param beta_vector Numeric vector containing values associated with the mesh vertices (e.g., posterior mean).
#' @param main_title Title of the plot.
#' @param resolution Resolution of the interpolation grid (default 100x100).
#' @return None. Displays a continuous thermal image plot.
plot_spde_field <- function(
  mesh,
  beta_vector,
  main_title = "Spatial Field",
  resolution = 100
) {
  # Extract spatial boundaries directly from the mesh
  xmin <- min(mesh$loc[, 1])
  xmax <- max(mesh$loc[, 1])
  ymin <- min(mesh$loc[, 2])
  ymax <- max(mesh$loc[, 2])

  # Create a dense grid for continuous spatial interpolation
  grid_coords <- expand.grid(
    x = seq(xmin, xmax, length.out = resolution),
    y = seq(ymin, ymax, length.out = resolution)
  )

  # Project vertex values onto the continuous grid using 1st-order FEM basis functions
  proj_vals <- fm_evaluate(
    mesh,
    loc = as.matrix(grid_coords),
    field = beta_vector
  )
  z_matrix <- matrix(proj_vals, nrow = resolution, ncol = resolution)

  # Setup color palette and render image
  col_palette <- hcl.colors(100, palette = "Spectral", rev = TRUE)

  image(
    seq(xmin, xmax, length.out = resolution),
    seq(ymin, ymax, length.out = resolution),
    z_matrix,
    col = col_palette,
    asp = 1,
    main = main_title,
    xlab = "X",
    ylab = "Y"
  )
}
