#' Load and Clean Annual Climate Data
#'
#' Loads the satellite temperature dataset (`climate_data.RData`), extracts
#' annual mean temperatures (`temp_tot`), and resolves geographic duplicates
#' at the poles (`latitude = +-90`) and along the `+-180` longitude seam,
#' yielding the N = 64,442 unique physical locations per year.
#'
#' @param path Character string. Path to `climate_data.RData`.
#' @param years Optional integer vector of years to keep (default `NULL` keeps all).
#' @param temp_col Character string. Name of the temperature column (default `"temp_tot"`).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{locations}: Data frame (`N x 2`) with columns `latitude`, `longitude`.
#'   \item \code{y_mat}: Numeric matrix (`N x M`) of annual temperatures (rows = locations, cols = years).
#'   \item \code{years}: Integer vector of length `M` with the included years.
#' }
#' @export
load_clean_climate_data <- function(
    path = "data/climate_data.RData",
    years = NULL,
    temp_col = "temp_tot"
) {
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }

  # Load into isolated environment to avoid polluting .GlobalEnv
  env <- new.env(parent = emptyenv())
  obj_names <- load(path, envir = env)

  if ("data" %in% obj_names && is.data.frame(env$data)) {
    df <- env$data
  } else {
    df_candidates <- Filter(function(nm) is.data.frame(env[[nm]]), obj_names)
    if (length(df_candidates) == 0) {
      stop("No data.frame found inside: ", path)
    }
    df <- env[[df_candidates[1]]]
  }

  req_cols <- c("latitude", "longitude", "year", temp_col)
  missing_cols <- setdiff(req_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  df <- df[, req_cols, drop = FALSE]
  names(df)[4] <- "temp"

  if (!is.null(years)) {
    df <- df[df$year %in% years, , drop = FALSE]
    if (nrow(df) == 0) stop("No observations left after filtering 'years'.")
  }

  if (anyNA(df$temp)) {
    warning("NA values detected in temperature column; removing incomplete rows.")
    df <- df[!is.na(df$temp), , drop = FALSE]
  }

  # Fast deduplication of physical locations on S^2:
  # 1) Interior points (|lat| < 90, lon not in {-180, 180}): already unique
  # 2) Seam points (|lat| < 90, lon in {-180, 180}): merge +180 into -180 by mean
  # 3) Poles (|lat| == 90): merge all longitudes into lon = 0 by mean
  is_pole <- abs(df$latitude) == 90L
  is_seam <- !is_pole & (df$longitude == 180L | df$longitude == -180L)
  is_interior <- !is_pole & !is_seam

  parts <- list()
  if (any(is_interior)) {
    parts$interior <- df[is_interior, , drop = FALSE]
  }
  if (any(is_seam)) {
    seam_df <- df[is_seam, , drop = FALSE]
    seam_df$longitude <- -180L
    parts$seam <- stats::aggregate(
      temp ~ latitude + longitude + year,
      data = seam_df,
      FUN = mean
    )
  }
  if (any(is_pole)) {
    pole_df <- df[is_pole, , drop = FALSE]
    pole_df$longitude <- 0L
    parts$poles <- stats::aggregate(
      temp ~ latitude + longitude + year,
      data = pole_df,
      FUN = mean
    )
  }

  clean_df <- do.call(rbind, parts)
  rownames(clean_df) <- NULL

  # Order deterministically: by year, then latitude, then longitude
  ord <- order(clean_df$year, clean_df$latitude, clean_df$longitude)
  clean_df <- clean_df[ord, , drop = FALSE]

  years_vec <- sort(unique(clean_df$year))
  M <- length(years_vec)

  # Build unique location table and N x M observation matrix
  loc_key <- paste(clean_df$latitude, clean_df$longitude, sep = "_")
  first_year_idx <- clean_df$year == years_vec[1]
  locations <- clean_df[first_year_idx, c("latitude", "longitude"), drop = FALSE]
  rownames(locations) <- NULL
  N <- nrow(locations)

  if (nrow(clean_df) == N * M &&
      all(loc_key == rep(loc_key[first_year_idx], times = M))) {
    # Fast vectorized reshape for complete panel
    y_mat <- matrix(clean_df$temp, nrow = N, ncol = M)
  } else {
    # General safe fallback if any location-year pair is missing
    uniq_keys <- unique(loc_key)
    first_occ <- match(uniq_keys, loc_key)
    locations <- clean_df[first_occ, c("latitude", "longitude"), drop = FALSE]
    rownames(locations) <- NULL
    N <- nrow(locations)
    y_mat <- matrix(NA_real_, nrow = N, ncol = M)
    r_idx <- match(loc_key, uniq_keys)
    c_idx <- match(clean_df$year, years_vec)
    y_mat[cbind(r_idx, c_idx)] <- clean_df$temp
  }

  colnames(y_mat) <- as.character(years_vec)

  list(
    locations = locations,
    y_mat = y_mat,
    years = years_vec
  )
}


#' Compute Great-Circle Distance (Haversine) in Kilometers
#'
#' @param lat,lon Numeric vectors of coordinates in degrees.
#' @param lat0,lon0 Numeric scalars of the reference point in degrees.
#' @param R_earth Earth radius in km (default `6371.0088`).
#'
#' @return Numeric vector of geodesic distances in km.
#' @keywords internal
haversine_km <- function(lat, lon, lat0, lon0, R_earth = 6371.0088) {
  deg2rad <- pi / 180
  phi <- lat * deg2rad
  phi0 <- lat0 * deg2rad
  dphi <- (lat - lat0) * deg2rad
  dlam <- (lon - lon0) * deg2rad

  a <- sin(dphi / 2)^2 + cos(phi0) * cos(phi) * sin(dlam / 2)^2
  c_ang <- 2 * atan2(sqrt(a), sqrt(pmax(0, 1 - a)))
  R_earth * c_ang
}


#' Project Spherical (lat, lon) Coordinates onto a Local Tangent Plane
#'
#' Maps geographic coordinates in degrees to planar `(x, y)` coordinates in km
#' centered at `(lat0, lon0)` without external GIS dependencies.
#'
#' @param lat,lon Numeric vectors of latitude and longitude in degrees.
#' @param lat0,lon0 Numeric scalars defining the projection center in degrees.
#' @param projection Character string: `"aeqd"` (Azimuthal Equidistant, exact
#'   radial distances from center) or `"laea"` (Lambert Azimuthal Equal-Area).
#' @param R_earth Earth radius in km (default `6371.0088`).
#'
#' @return Numeric matrix (`N x 2`) with columns `x` (East-West, km) and `y` (North-South, km).
#' @export
project_to_plane <- function(
    lat,
    lon,
    lat0,
    lon0,
    projection = c("aeqd", "laea"),
    R_earth = 6371.0088
) {
  projection <- match.arg(projection)
  deg2rad <- pi / 180

  phi <- lat * deg2rad
  phi0 <- lat0 * deg2rad
  dlam <- (lon - lon0) * deg2rad

  cos_phi <- cos(phi)
  sin_phi <- sin(phi)
  cos_phi0 <- cos(phi0)
  sin_phi0 <- sin(phi0)
  cos_dlam <- cos(dlam)
  sin_dlam <- sin(dlam)

  if (projection == "aeqd") {
    a <- sin((phi - phi0) / 2)^2 + cos_phi0 * cos_phi * sin(dlam / 2)^2
    c_ang <- 2 * atan2(sqrt(a), sqrt(pmax(0, 1 - a)))
    k <- rep(1.0, length(c_ang))
    nz <- c_ang >= 1e-12 & c_ang < (pi - 1e-9)
    k[nz] <- c_ang[nz] / sin(c_ang[nz])
    k[c_ang >= (pi - 1e-9)] <- NA_real_
  } else {
    # Lambert Azimuthal Equal-Area
    denom <- 1 + sin_phi0 * sin_phi + cos_phi0 * cos_phi * cos_dlam
    k <- ifelse(denom > 1e-12, sqrt(2 / denom), NA_real_)
  }

  x <- R_earth * k * cos_phi * sin_dlam
  y <- R_earth * k * (cos_phi0 * sin_phi - sin_phi0 * cos_phi * cos_dlam)

  cbind(x = x, y = y)
}


#' Extract a Local Planar Region from Cleaned Climate Data
#'
#' Subsets the global climate dataset around a reference coordinate using either
#' a planar geometric shape (`"square"`, `"rectangle"`, `"diamond"`, `"circle"`)
#' or geodesic K-Nearest Neighbors (`"knn"`), and returns flattened coordinates,
#' observations, and domain boundary vertices.
#'
#' @param clean_data Output list from [load_clean_climate_data()].
#' @param center Numeric vector `c(lat, lon)` in degrees for the center of the region.
#' @param origin Optional numeric vector `c(lat, lon)` in degrees specifying the
#'   South-West corner (for `"square"` and `"rectangle"` only). If provided,
#'   overrides `center` and defines the box `[0, width_km] x [0, height_km]`.
#' @param shape Character string. One of `"square"`, `"rectangle"`, `"diamond"`,
#'   `"circle"`, or `"knn"`.
#' @param side_km Numeric scalar. Side length in km (for `"square"`, or regular `"diamond"`).
#' @param width_km,height_km Numeric scalars in km (for `"rectangle"`, or diagonals for `"diamond"`).
#' @param radius_km Numeric scalar. Radius in km (for `"circle"`).
#' @param k Integer. Number of nearest neighbors (for `"knn"`).
#' @param projection Character string. `"aeqd"` (default) or `"laea"`.
#' @param scale_to_unit Logical. If `TRUE`, isotropically rescales `(x, y)` and
#'   `boundary` by the characteristic domain length (preserving anisotropy) so
#'   coordinates live on a ~unit scale compatible with unit-square PC priors.
#' @param n_circle_pts Integer. Number of vertices used to represent a circular boundary.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{y_mat}: Matrix (`N_sub x M`) of annual temperatures.
#'   \item \code{coords}: Matrix (`N_sub x 2`) of planar coordinates `(x, y)`.
#'   \item \code{latlon}: Data frame (`N_sub x 2`) of original `(latitude, longitude)`.
#'   \item \code{boundary}: Matrix (`B x 2`) of CCW-ordered boundary vertices.
#'   \item \code{years}: Integer vector of years (`length M`).
#'   \item \code{meta}: List of metadata (shape, center/origin, projection, units, scale_factor).
#' }
#' @export
extract_planar_region <- function(
    clean_data,
    center = c(41.9028, 12.4964),
    origin = NULL,
    shape = c("square", "rectangle", "diamond", "circle", "knn"),
    side_km = 2000,
    width_km = NULL,
    height_km = NULL,
    radius_km = 1000,
    k = 400L,
    projection = c("aeqd", "laea"),
    scale_to_unit = FALSE,
    n_circle_pts = 64L
) {
  shape <- match.arg(shape)
  projection <- match.arg(projection)

  locs <- clean_data$locations
  use_origin <- !is.null(origin)

  if (use_origin) {
    if (!(shape %in% c("square", "rectangle"))) {
      stop("'origin' (South-West corner) is only supported for 'square' and 'rectangle'.")
    }
    ref_pt <- origin
  } else {
    ref_pt <- center
  }

  if (length(ref_pt) != 2 || !is.numeric(ref_pt)) {
    stop("Reference point ('center' or 'origin') must be a numeric vector c(lat, lon).")
  }
  lat0 <- ref_pt[1]
  lon0 <- ref_pt[2]

  # Project all locations onto local tangent plane (in km)
  xy_all <- project_to_plane(
    lat = locs$latitude,
    lon = locs$longitude,
    lat0 = lat0,
    lon0 = lon0,
    projection = projection
  )
  valid <- stats::complete.cases(xy_all)

  if (shape == "square") {
    W <- side_km
    H <- side_km
    if (use_origin) {
      keep <- valid & xy_all[, 1] >= 0 & xy_all[, 1] <= W &
        xy_all[, 2] >= 0 & xy_all[, 2] <= H
      boundary <- rbind(c(0, 0), c(W, 0), c(W, H), c(0, H))
    } else {
      keep <- valid & abs(xy_all[, 1]) <= W / 2 & abs(xy_all[, 2]) <= H / 2
      boundary <- rbind(
        c(-W / 2, -H / 2),
        c( W / 2, -H / 2),
        c( W / 2,  H / 2),
        c(-W / 2,  H / 2)
      )
    }
    char_scale <- side_km

  } else if (shape == "rectangle") {
    W <- if (is.null(width_km)) side_km else width_km
    H <- if (is.null(height_km)) side_km / 2 else height_km
    if (use_origin) {
      keep <- valid & xy_all[, 1] >= 0 & xy_all[, 1] <= W &
        xy_all[, 2] >= 0 & xy_all[, 2] <= H
      boundary <- rbind(c(0, 0), c(W, 0), c(W, H), c(0, H))
    } else {
      keep <- valid & abs(xy_all[, 1]) <= W / 2 & abs(xy_all[, 2]) <= H / 2
      boundary <- rbind(
        c(-W / 2, -H / 2),
        c( W / 2, -H / 2),
        c( W / 2,  H / 2),
        c(-W / 2,  H / 2)
      )
    }
    char_scale <- max(W, H)

  } else if (shape == "diamond") {
    # Half-diagonals a (East-West) and b (North-South)
    a <- if (is.null(width_km)) side_km / sqrt(2) else width_km / 2
    b <- if (is.null(height_km)) side_km / sqrt(2) else height_km / 2
    keep <- valid & (abs(xy_all[, 1]) / a + abs(xy_all[, 2]) / b <= 1 + 1e-12)
    boundary <- rbind(c(a, 0), c(0, b), c(-a, 0), c(0, -b))
    char_scale <- max(2 * a, 2 * b)

  } else if (shape == "circle") {
    # Use exact geodesic distance so circle matches spherical cap
    dist_km <- haversine_km(locs$latitude, locs$longitude, lat0, lon0)
    keep <- valid & (dist_km <= radius_km)
    theta <- seq(0, 2 * pi, length.out = n_circle_pts + 1L)[- (n_circle_pts + 1L)]
    boundary <- cbind(x = radius_km * cos(theta), y = radius_km * sin(theta))
    char_scale <- 2 * radius_km

  } else if (shape == "knn") {
    k <- as.integer(k)
    if (k < 3L || k > nrow(locs)) {
      stop("'k' must be between 3 and ", nrow(locs), ".")
    }
    dist_km <- haversine_km(locs$latitude, locs$longitude, lat0, lon0)
    ord_d <- order(dist_km)
    keep_idx <- ord_d[seq_len(k)]
    keep <- rep(FALSE, nrow(locs))
    keep[keep_idx] <- TRUE

    pts_sub <- xy_all[keep_idx, , drop = FALSE]
    # Convex hull in counter-clockwise order (grDevices::chull returns clockwise)
    hull_idx <- rev(grDevices::chull(pts_sub[, 1], pts_sub[, 2]))
    boundary <- pts_sub[hull_idx, , drop = FALSE]
    char_scale <- 2 * max(dist_km[keep_idx])
  }

  idx <- which(keep)
  if (length(idx) < 3L) {
    stop(
      "Too few grid points (", length(idx), ") inside the requested region. ",
      "Remember the grid resolution is 1 deg (~111 km): try increasing the region size."
    )
  }

  coords <- xy_all[idx, , drop = FALSE]
  colnames(boundary) <- c("x", "y")
  latlon <- locs[idx, , drop = FALSE]
  rownames(latlon) <- NULL
  y_mat <- clean_data$y_mat[idx, , drop = FALSE]

  # Optional isotropic scaling (divides both x and y by the SAME factor)
  scale_factor <- 1.0
  units_label <- "km"
  if (scale_to_unit) {
    scale_factor <- char_scale
    coords <- coords / scale_factor
    boundary <- boundary / scale_factor
    units_label <- "unit_scaled"
  }

  list(
    y_mat = y_mat,
    coords = coords,
    latlon = latlon,
    boundary = boundary,
    years = clean_data$years,
    meta = list(
      shape = shape,
      ref_type = if (use_origin) "origin_sw" else "center",
      ref_latlon = c(lat = lat0, lon = lon0),
      projection = projection,
      units = units_label,
      scale_factor_km = scale_factor,
      n_obs = length(idx),
      n_years = ncol(y_mat)
    )
  )
}


#' Build Planar Triangular Mesh and FEM Evaluation Matrix for a Region
#'
#' Creates an `fm_mesh_2d` object (via `fmesher`) covering the extracted planar
#' region and computes the sparse basis evaluation matrix `Phi` (`N_obs x V`)
#' mapping mesh nodes to observation locations.
#'
#' @param region Output list from [extract_planar_region()].
#' @param max_edge Numeric scalar or length-2 vector `c(inner, outer)` controlling
#'   maximum triangle edge length (in the same units as `region$coords`). If `NULL`,
#'   a sensible default based on domain span is chosen.
#' @param offset Numeric scalar (outer extension width for SPDE boundary-effect
#'   mitigation, e.g. 20% of domain span). Set to `NULL` or `0` for no outer extension.
#' @param cutoff Minimum allowed distance between mesh nodes. If `NULL`, set to
#'   `0.2 * max_edge[1]`.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{mesh}: An `fm_mesh_2d` object compatible with [build_spde_matrices()].
#'   \item \code{Phi}: Sparse matrix (`N_obs x V`, class `dgCMatrix`) evaluating
#'     degree-1 FEM basis functions at `region$coords`.
#' }
#' @export
build_planar_mesh <- function(
    region,
    max_edge = NULL,
    offset = NULL,
    cutoff = NULL
) {
  if (!requireNamespace("fmesher", quietly = TRUE)) {
    stop("Package 'fmesher' is required to build the mesh.")
  }

  span <- max(
    diff(range(region$boundary[, 1])),
    diff(range(region$boundary[, 2]))
  )

  use_extension <- !is.null(offset) && any(offset > 0)

  if (is.null(max_edge)) {
    max_edge <- if (use_extension) c(span / 12, span / 4) else span / 12
  }
  if (is.null(cutoff)) {
    cutoff <- 0.2 * max_edge[1]
  }

  # All supported shapes (square, rectangle, diamond, circle, knn convex hull)
  # are convex polygons, so loc.domain = region$boundary builds the exact domain cleanly.
  mesh <- if (use_extension) {
    fmesher::fm_mesh_2d_inla(
      loc.domain = region$boundary,
      max.edge = max_edge,
      offset = offset,
      cutoff = cutoff
    )
  } else {
    fmesher::fm_mesh_2d_inla(
      loc.domain = region$boundary,
      max.edge = max_edge[1],
      cutoff = cutoff
    )
  }

  # Sparse N_obs x V evaluation matrix Phi
  Phi <- fmesher::fm_basis(mesh, loc = region$coords)

  list(
    mesh = mesh,
    Phi = Phi
  )
}


#' End-to-End Preparation of Real Planar Climate Data for SPDE Testing
#'
#' Convenience wrapper combining data loading/cleaning, local planar extraction,
#' and optional mesh + `Phi` matrix construction.
#'
#' @param path Path to `climate_data.RData`.
#' @param years Optional integer vector of years to include.
#' @param build_mesh Logical. If `TRUE` (default), appends `mesh` and `Phi`.
#' @param max_edge,offset,cutoff Mesh parameters passed to [build_planar_mesh()].
#' @param ... Arguments passed to [extract_planar_region()] (e.g., `center`,
#'   `origin`, `shape`, `side_km`, `k`, `scale_to_unit`).
#'
#' @return A list with `y_mat`, `coords`, `latlon`, `boundary`, `years`, `meta`,
#'   and (if `build_mesh = TRUE`) `mesh` and `Phi`.
#' @export
prepare_real_planar_data <- function(
    path = "data/climate_data.RData",
    years = NULL,
    build_mesh = TRUE,
    max_edge = NULL,
    offset = NULL,
    cutoff = NULL,
    ...
) {
  clean_data <- load_clean_climate_data(path = path, years = years)
  res <- extract_planar_region(clean_data = clean_data, ...)

  if (build_mesh) {
    mesh_out <- build_planar_mesh(
      region = res,
      max_edge = max_edge,
      offset = offset,
      cutoff = cutoff
    )
    res$mesh <- mesh_out$mesh
    res$Phi <- mesh_out$Phi
  }

  res
}
