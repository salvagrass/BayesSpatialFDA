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


#' Anisotropy matrix entries H_v from the half-angle vector v.
#'
#' Implements H_v = cosh(||v||) I + sinh(||v||)/||v|| * [[v1 v2],[v2 -v1]],
#' with the removable isotropic limit sinh(||v||)/||v|| -> 1 as ||v|| -> 0.
#'
#' @param v1,v2 Components of the anisotropy vector.
#' @return Named list with scalars H11, H22, H12.
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


#' Closed boundary loops of a mesh
#'
#' Chains the boundary edges of the mesh (the triangle edges that belong to
#' exactly one triangle) into closed loops of vertex indices: the outer
#' boundary and one loop per hole.
#'
#' This function is used internally as a helper by `plot_spde_field()`.
#'
#' @param mesh The mesh object created with fmesher.
#' @return List of integer vectors of mesh vertex indices; each is one loop,
#'   without repeating its first vertex at the end.
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


#' Plot a continuous SPDE spatial field using mesh and 1st-order FEM basis functions
#'
#' @description Plots the piecewise-linear FEM field defined by the mesh vertex
#'   values, with its own colour legend. The field is evaluated exactly at the
#'   pixel centres and the picture is cut along the
#'   true mesh boundary. Base graphics: it draws into the current panel, so it
#'   can be combined with `par(mfrow = ...)` or `layout()`.
#' @param mesh The mesh object created with fmesher.
#' @param beta_vector Numeric vector containing values associated with the mesh vertices (e.g., posterior mean).
#' @param main_title Title of the plot (a string or a plotmath expression).
#' @param resolution Resolution of the interpolation grid (default 100x100).
#' @return Invisibly, a list with the pixel centres `x`, `y` and the matrix `z`
#'   of plotted values (`NA` outside the mesh). Called for its side effect.
plot_spde_field <- function(
  mesh,
  beta_vector,
  main_title = "Spatial Field",
  resolution = 100
) {
  beta_vector <- as.numeric(beta_vector)
  loc <- mesh$loc[, 1:2, drop = FALSE]
  xr <- range(loc[, 1])
  yr <- range(loc[, 2])
  w <- diff(xr)
  h <- diff(yr)

  # Pixel centres: resolution x resolution cells tiling the mesh bounding box
  gx <- xr[1] + (seq_len(resolution) - 0.5) * w / resolution
  gy <- yr[1] + (seq_len(resolution) - 0.5) * h / resolution
  grid_coords <- as.matrix(expand.grid(x = gx, y = gy))

  # Evaluate the 1st-order FEM field at the pixel centres (NA outside the mesh)
  z <- matrix(
    fm_evaluate(mesh, loc = grid_coords, field = beta_vector),
    nrow = resolution,
    ncol = resolution
  )

  # Pixels straddling the boundary whose centre falls outside the mesh would be
  # left blank and give a jagged edge. They take the field value at the nearest
  # boundary point (exact: the field is linear along an edge); the part of them
  # lying outside the mesh is masked below.
  grow <- function(m) {
    padded <- matrix(FALSE, nrow(m) + 2, ncol(m) + 2)
    padded[1 + seq_len(nrow(m)), 1 + seq_len(ncol(m))] <- m
    out <- m
    for (di in 0:2) {
      for (dj in 0:2) {
        out <- out | padded[di + seq_len(nrow(m)), dj + seq_len(ncol(m))]
      }
    }
    out
  }
  loops <- mesh_boundary_loops(mesh)
  edge <- which(is.na(z) & grow(grow(!is.na(z))))
  if (length(edge) > 0) {
    seg <- do.call(
      rbind,
      lapply(loops, function(loop) cbind(loop, c(loop[-1], loop[1])))
    )
    a <- loc[seg[, 1], , drop = FALSE]
    d <- loc[seg[, 2], , drop = FALSE] - a
    px <- outer(grid_coords[edge, 1], a[, 1], "-") # pixel x segment
    py <- outer(grid_coords[edge, 2], a[, 2], "-")
    dx <- rep(d[, 1], each = length(edge))
    dy <- rep(d[, 2], each = length(edge))
    t <- pmin(pmax((px * dx + py * dy) / (dx^2 + dy^2), 0), 1)
    dist2 <- (px - t * dx)^2 + (py - t * dy)^2
    nearest <- max.col(-dist2, ties.method = "first")
    t <- t[cbind(seq_along(edge), nearest)]
    z[edge] <- (1 - t) *
      beta_vector[seg[nearest, 1]] +
      t * beta_vector[seg[nearest, 2]]
  }

  # A P1 field attains its extremes at the vertices
  zlim <- range(beta_vector, finite = TRUE)
  if (diff(zlim) == 0) {
    zlim <- zlim + c(-0.5, 0.5)
  }
  col_palette <- hcl.colors(256, palette = "Spectral", rev = TRUE)
  ink <- "grey25"

  # Legend geometry: a bar as tall as the map, to its right, inside the plot
  # region. The plot region is sized so that map + legend fit at aspect 1.
  ticks <- pretty(zlim, n = 5)
  ticks <- ticks[ticks >= zlim[1] & ticks <= zlim[2]]
  tick_labels <- format(ticks, trim = TRUE)
  cex_axis <- par("cex.axis")
  bar_gap <- 0.04 * w
  bar_width <- 0.045 * w

  plot.new()
  inch <- par("cin")[2] * par("cex") # one text line, in inches
  label_inches <- max(strwidth(tick_labels, units = "inches", cex = cex_axis)) +
    0.6 * inch
  pin <- par("pin")
  scale <- min(
    (pin[1] - label_inches) / (w + bar_gap + bar_width),
    pin[2] / h
  ) # inches per user unit
  plot.window(
    xlim = c(xr[1], xr[2] + bar_gap + bar_width + label_inches / scale),
    ylim = yr,
    asp = 1,
    xaxs = "i",
    yaxs = "i"
  )
  line <- inch / scale # one text line, in user units

  image(gx, gy, z, zlim = zlim, col = col_palette, add = TRUE, useRaster = TRUE)

  # Mask everything outside the mesh (even-odd rule: bounding box minus the
  # mesh domain, holes included), then outline the boundary.
  pad <- 0.01 * max(w, h)
  mask_x <- c(xr[1] - pad, xr[2] + pad, xr[2] + pad, xr[1] - pad)
  mask_y <- c(yr[1] - pad, yr[1] - pad, yr[2] + pad, yr[2] + pad)
  for (loop in loops) {
    mask_x <- c(mask_x, NA, loc[loop, 1])
    mask_y <- c(mask_y, NA, loc[loop, 2])
  }
  bg <- par("bg")
  if (grDevices::col2rgb(bg, alpha = TRUE)[4] == 0) {
    bg <- "white"
  }
  polypath(mask_x, mask_y, rule = "evenodd", col = bg, border = NA)
  for (loop in loops) {
    polygon(loc[loop, 1], loc[loop, 2], border = ink, lwd = 0.8)
  }

  # Axes attached to the map (ticks and numbers only)
  x_at <- pretty(xr)
  y_at <- pretty(yr)
  axis(
    1,
    at = x_at[x_at >= xr[1] & x_at <= xr[2]],
    pos = yr[1] - 0.3 * line,
    lwd = 0,
    lwd.ticks = 0.8,
    tcl = -0.25,
    mgp = c(3, 0.3, 0),
    col.ticks = ink,
    col.axis = ink
  )
  axis(
    2,
    at = y_at[y_at >= yr[1] & y_at <= yr[2]],
    pos = xr[1] - 0.3 * line,
    lwd = 0,
    lwd.ticks = 0.8,
    tcl = -0.25,
    mgp = c(3, 0.5, 0),
    las = 1,
    col.ticks = ink,
    col.axis = ink
  )

  # Colour legend
  bar_x <- xr[2] + bar_gap + c(0, bar_width)
  rasterImage(
    as.raster(matrix(rev(col_palette), ncol = 1)),
    bar_x[1],
    yr[1],
    bar_x[2],
    yr[2],
    interpolate = FALSE
  )
  rect(bar_x[1], yr[1], bar_x[2], yr[2], border = ink, lwd = 0.8)
  tick_y <- yr[1] + (ticks - zlim[1]) / diff(zlim) * h
  tick_x <- bar_x[2] + 0.25 * line
  segments(bar_x[2], tick_y, tick_x, tick_y, col = ink, lwd = 0.8)
  text(
    bar_x[2] + 0.5 * line,
    tick_y,
    tick_labels,
    adj = c(0, 0.5),
    cex = cex_axis,
    col = ink,
    xpd = NA
  )

  # Title centred on the map, just above it
  text(
    mean(xr),
    yr[2] + 0.7 * line,
    main_title,
    adj = c(0.5, 0),
    cex = par("cex.main"),
    font = par("font.main"),
    col = par("col.main"),
    xpd = NA
  )

  invisible(list(x = gx, y = gy, z = z))
}

