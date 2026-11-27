#' Sample from a Multivariate Normal Distribution
#'
#' Generates a single random sample from a multivariate normal distribution.
#' It handles scalar inputs for both mean and covariance/precision efficiently
#' using sparse matrices.
#'
#' @param mean A numeric vector representing the mean. If a scalar is provided,
#'   it will be repeated `N` times.
#' @param mat A numeric or sparse matrix representing the covariance matrix or
#'   its inverse (precision matrix). If a scalar is provided, it is treated as
#'   the variance (or precision) of independent variables and an internal sparse
#'   diagonal matrix is generated.
#' @param N An optional integer specifying the dimension of the vector. If `NULL`,
#'   it is inferred from the length of `mean` or the number of rows of `mat`.
#' @param is_precision Logical. If `TRUE`, `mat` is treated as a precision matrix
#'   rather than a covariance matrix. Default is `FALSE`.
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


#' Anisotropy Matrix Entries from the Half-Angle Vector
#'
#' Computes the entries of the symmetric positive-definite diffusion tensor
#' \eqn{H_v = \cosh(\|v\|) I_2 + \frac{\sinh(\|v\|)}{\|v\|} \begin{pmatrix} v_1 & v_2 \\ v_2 & -v_1 \end{pmatrix}},
#' using the removable isotropic limit \eqn{\sinh(\|v\|)/\|v\| \to 1} as
#' \eqn{\|v\| \to 0}.
#'
#' @param v1 Numeric scalar. First component of the anisotropy vector.
#' @param v2 Numeric scalar. Second component of the anisotropy vector.
#'
#' @return A named list with scalars `H11`, `H22`, and `H12`.
anisotropy_from_v <- function(v1, v2) {
  norm_v <- sqrt(v1^2 + v2^2)
  if (norm_v < .Machine$double.eps) {
    # Isotropic limit: sinh(||v||)/||v|| -> 1 as ||v|| -> 0, H_v -> I
    sinhc_v <- 1.0
  } else {
    sinhc_v <- sinh(norm_v) / norm_v
  }

  list(
    H11 = cosh(norm_v) + sinhc_v * v1,
    H22 = cosh(norm_v) - sinhc_v * v1,
    H12 = sinhc_v * v2
  )
}


#' Closed Boundary Loops of a Mesh
#'
#' Chains the boundary edges of the mesh (the triangle edges that belong to
#' exactly one triangle) into closed loops of vertex indices: the outer
#' boundary and one loop per hole.
#'
#' This function is used internally as a helper by [plot_spde_field()].
#'
#' @param mesh A mesh object created with `fmesher`.
#'
#' @return A list of integer vectors of mesh vertex indices; each element is
#'   one closed loop without repeating its first vertex at the end.
mesh_boundary_loops <- function(mesh) {
  tv <- mesh$graph$tv
  edges <- rbind(tv[, c(1, 2)], tv[, c(2, 3)], tv[, c(3, 1)])
  key <- paste(pmin(edges[, 1], edges[, 2]), pmax(edges[, 1], edges[, 2]))
  seg <- edges[key %in% names(which(table(key) == 1)), , drop = FALSE]

  used <- logical(nrow(seg))
  loops <- list()
  while (!all(used)) {
    i <- which(!used)[1]
    loop <- integer(0)
    while (length(i) == 1) {
      used[i] <- TRUE
      loop <- c(loop, seg[i, 1])
      i <- which(!used & seg[, 1] == seg[i, 2])[1]
      if (is.na(i)) i <- integer(0)
    }
    loops[[length(loops) + 1]] <- loop
  }
  loops
}


#' Plot a Continuous SPDE Spatial Field
#'
#' @description Plots the piecewise-linear FEM field defined by the mesh vertex
#'   values, with its own colour legend. The field is evaluated at the pixel
#'   centres and clipped along the true mesh boundary using base graphics.
#'
#' @param mesh A mesh object created with `fmesher`.
#' @param beta_vector Numeric vector containing values associated with the mesh
#'   vertices (e.g., posterior mean).
#' @param main_title Character string or `plotmath` expression for the plot title.
#'   Default is `"Spatial Field"`.
#' @param resolution Integer. Resolution of the interpolation grid
#'   (`resolution x resolution`). Default is `100`.
#'
#' @return Invisibly, a list with the pixel centres `x`, `y` and the matrix `z`
#'   of plotted values (`NA` outside the mesh). Called for its side effect.
plot_spde_field <- function(
  mesh,
  beta_vector,
  main_title = "Spatial Field",
  resolution = 100,
  xlim = c(0, 10),
  ylim = c(0, 10),
  zlim = NULL
) {
  grid_coords <- expand.grid(
    x = seq(xlim[1], xlim[2], length.out = resolution),
    y = seq(ylim[1], ylim[2], length.out = resolution)
  )

  proj_vals <- fm_evaluate(
    mesh,
    loc = as.matrix(grid_coords),
    field = beta_vector
  )
  z_matrix <- matrix(proj_vals, nrow = resolution, ncol = resolution)

  col_palette <- hcl.colors(100, palette = "Spectral", rev = TRUE)

  if (is.null(zlim)) {
    zlim <- range(z_matrix, na.rm = TRUE)
  }

  image(
    seq(xlim[1], xlim[2], length.out = resolution),
    seq(ylim[1], ylim[2], length.out = resolution),
    z_matrix,
    col = col_palette,
    zlim = zlim,
    asp = 1,
    main = main_title,
    xlab = "X",
    ylab = "Y"
  )
}
#' Build the Anisotropic SPDE Precision Matrix
#'
#' @description
#' Assembles the sparse Gaussian Markov Random Field (GMRF) precision matrix
#' \eqn{Q(\theta)} for the anisotropic Matérn SPDE with smoothness \eqn{\nu = 1}
#' under the half-angle parameterization:
#' \deqn{Q(\theta) = \frac{1}{4\pi\kappa^2\sigma_\beta^2} (\kappa^2 \tilde{C} + G_{H_v}) \tilde{C}^{-1} (\kappa^2 \tilde{C} + G_{H_v})}
#' where \eqn{\tilde{C}} is the lumped diagonal mass matrix and \eqn{G_{H_v}}
#' is the anisotropic stiffness matrix built from the precomputed geometric
#' stiffness components and [anisotropy_from_v()].
#'
#' @param spde_matrices A named list of precomputed finite element matrices
#'   containing the mass matrix `C` and the stiffness components `K11`, `K12`,
#'   and `K22`.
#' @param kappa Numeric scalar (\eqn{\kappa > 0}). SPDE spatial inverse
#'   length-scale parameter, related to the empirical correlation range via
#'   \eqn{\rho = \sqrt{8}/\kappa}.
#' @param sigma_beta Numeric scalar (\eqn{\sigma_\beta > 0}). Marginal
#'   standard deviation of the latent spatial field.
#' @param v1 Numeric scalar. First component of the unconstrained anisotropy
#'   vector \eqn{v = (v_1, v_2)^\top}. Default is `0.0`.
#' @param v2 Numeric scalar. Second component of the unconstrained anisotropy
#'   vector \eqn{v = (v_1, v_2)^\top}. Default is `0.0`.
#'
#' @return A symmetric sparse precision matrix of dimension \eqn{V \times V}.
#' @importFrom Matrix Diagonal rowSums
#' @export
build_precision_matrix <- function(
  spde_matrices,
  kappa,
  sigma_beta,
  v1 = 0.0,
  v2 = 0.0
) {
  H <- anisotropy_from_v(v1, v2)

  global_K <- H$H11 *
    spde_matrices$K11 +
    H$H12 * spde_matrices$K12 +
    H$H22 * spde_matrices$K22

  c_lumped <- rowSums(spde_matrices$C)
  C_tilde <- Diagonal(x = c_lumped)
  C_tilde_inv <- Diagonal(x = 1 / c_lumped)

  core_mat <- global_K + (kappa^2) * C_tilde
  scale <- 1 / (4 * pi * kappa^2 * sigma_beta^2)

  scale * (core_mat %*% C_tilde_inv %*% core_mat)
}
