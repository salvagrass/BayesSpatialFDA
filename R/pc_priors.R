# R/pc_priors.R
# Penalized-complexity (PC) priors for the anisotropic SPDE hyperparameters
# theta = (v1, v2, kappa, sigma_beta): the distance function f, the calibration
# of the prior rates from tail quantiles, and the map Psi that pushes a standard
# Gaussian N(0, I_4) forward to the full PC prior.

# f(0) = sqrt(4 * pi / 3), the value of the PC distance at the isotropic model.
.pc_f0 <- sqrt(4 * pi / 3)

# acosh(1 + delta) for delta >= 0, written so that it stays accurate as
# delta -> 0 (acosh(1 + delta) itself cancels catastrophically there).
.acosh1p <- function(delta) {
  log1p(delta + sqrt(delta * (delta + 2)))
}

# Principal Lambert W evaluated at exp(L): the unique w > 0 with w + log(w) = L.
# Working with the log of the argument keeps the calibration of lambda_theta
# finite when exp(lambda_v * f(0)) would overflow (a0 close to 1).
# Newton's method on the concave increasing map w -> w + log(w) - L; the
# starting point satisfies log(w) < L, so every iterate stays positive.
.lambert_w0_exp <- function(L) {
  if (!is.finite(L)) {
    stop("'L' must be finite.")
  }
  if (L < -35) {
    # W0(x) = x - x^2 + O(x^3): the correction is below double precision here.
    return(exp(L))
  }
  w <- if (L > 1) L - log(L) else exp(L) / (1 + exp(L))
  for (iter in 1:100) {
    w_new <- w * (1 + L - log(w)) / (1 + w)
    # Convergence is quadratic (relative error after a step is below half the
    # squared relative step), so a step this small leaves only rounding error.
    converged <- abs(w_new - w) <= 1e-9 * w_new
    w <- w_new
    if (converged) {
      return(w)
    }
  }
  stop("Lambert W iteration did not converge.")
}

#' PC-prior distance function for the anisotropy vector
#'
#' `pc_f(r)` is the distance `f(r) = sqrt((pi / 3) * (3 * cosh(2 * r) + 1))`
#' appearing in the joint PC prior on `(kappa, v)`, evaluated at `r = ||v||`.
#' `pc_f_inv()` is its inverse on `[f(0), Inf)` and `pc_f_prime()` its
#' derivative `sqrt(pi) * sinh(2 * r) / sqrt(cosh(2 * r) + 1 / 3)`.
#'
#' The inverse is `r = acosh(x^2 / pi - 1 / 3) / 2`. It is evaluated as
#' `acosh(1 + delta)` with `delta = (x - f(0)) * (x + f(0)) / pi` through
#' `log1p`, which is the same function but does not lose accuracy for `x`
#' close to `f(0)`.
#'
#' @param r Numeric vector of non-negative anisotropy norms `||v||`.
#' @param x Numeric vector of distances, `x >= f(0) = sqrt(4 * pi / 3)`.
#'   Values below `f(0)` are outside the domain and give `NaN`.
#' @return A numeric vector of the same length as the input.
#' @export
pc_f <- function(r) {
  sqrt((pi / 3) * (3 * cosh(2 * r) + 1))
}

#' @rdname pc_f
#' @export
pc_f_inv <- function(x) {
  delta <- (x - .pc_f0) * (x + .pc_f0) / pi
  0.5 * .acosh1p(delta)
}

#' @rdname pc_f
#' @export
pc_f_prime <- function(r) {
  sqrt(pi) * sinh(2 * r) / sqrt(cosh(2 * r) + 1 / 3)
}

#' Calibrate the PC-prior rates from tail quantiles
#'
#' The rates of the PC prior on `theta = (v1, v2, kappa, sigma_beta)` are fixed
#' through three tail statements:
#'
#' * `P(exp(||v||) > a0) = alpha_v` (anisotropy ratio), giving
#'   `lambda_v = -log(alpha_v) / (f(log(a0)) - f(0))`;
#' * `P(kappa > kappa0) = alpha_rho`, i.e. `P(rho < rho0) = alpha_rho` for the
#'   correlation range `rho = sqrt(8) / kappa` with `kappa0 = sqrt(8) / rho0`,
#'   giving
#'   `lambda_theta = (W0(exp(lambda_v * f(0)) * lambda_v * f(0) / alpha_rho) / f(0) - lambda_v) / kappa0`
#'   with `W0` the principal branch of the Lambert W function;
#' * `P(sigma_beta > sigma0) = alpha_sigma`, giving
#'   `lambda_sigma_beta = -log(alpha_sigma) / sigma0`.
#'
#' The statement on `kappa` is marginal over `v`: conditionally on `v`,
#' `kappa` is exponential with rate `lambda_theta * f(||v||)`.
#'
#' @param a0 Surprisingly large anisotropy ratio, `a0 > 1`.
#' @param alpha_v Tail probability attached to `a0`, in `(0, 1)`.
#' @param kappa0 Surprisingly large `kappa` (`sqrt(8) / rho0` for a
#'   surprisingly small range `rho0`), `kappa0 > 0`.
#' @param lambda_v Anisotropy rate, as returned by `pc_calibrate_v()`.
#' @param alpha_rho Tail probability attached to `kappa0`, in `(0, 1)`.
#' @param sigma0 Surprisingly large marginal standard deviation, `sigma0 > 0`.
#' @param alpha_sigma Tail probability attached to `sigma0`, in `(0, 1)`.
#' @return `pc_calibrate_v()`, `pc_calibrate_theta()` and
#'   `pc_calibrate_sigma()` return the corresponding rate as a scalar.
#'   `pc_calibrate_all()` returns the list
#'   `list(lambda_v, lambda_theta, lambda_sigma_beta)` expected by
#'   `psi_reparam()`.
#' @name pc_calibrate
NULL

#' @rdname pc_calibrate
#' @export
pc_calibrate_v <- function(a0, alpha_v) {
  stopifnot(a0 > 1, alpha_v > 0, alpha_v < 1)
  -log(alpha_v) / (pc_f(log(a0)) - .pc_f0)
}

#' @rdname pc_calibrate
#' @export
pc_calibrate_theta <- function(kappa0, lambda_v, alpha_rho) {
  stopifnot(kappa0 > 0, lambda_v > 0, alpha_rho > 0, alpha_rho < 1)
  x <- lambda_v * .pc_f0
  # W0(exp(x) * x / alpha_rho), with the argument passed on the log scale
  w <- .lambert_w0_exp(x + log(x) - log(alpha_rho))
  (w / .pc_f0 - lambda_v) / kappa0
}

#' @rdname pc_calibrate
#' @export
pc_calibrate_sigma <- function(sigma0, alpha_sigma) {
  stopifnot(sigma0 > 0, alpha_sigma > 0, alpha_sigma < 1)
  -log(alpha_sigma) / sigma0
}

#' @rdname pc_calibrate
#' @export
pc_calibrate_all <- function(
  a0,
  alpha_v,
  kappa0,
  alpha_rho,
  sigma0,
  alpha_sigma
) {
  lambda_v <- pc_calibrate_v(a0, alpha_v)
  list(
    lambda_v = lambda_v,
    lambda_theta = pc_calibrate_theta(kappa0, lambda_v, alpha_rho),
    lambda_sigma_beta = pc_calibrate_sigma(sigma0, alpha_sigma)
  )
}

#' Map a standard Gaussian vector to the PC prior on the SPDE hyperparameters
#'
#' For `Y ~ N(0, I_4)` the image `Psi(Y) = (v1, v2, kappa, sigma_beta)` is
#' distributed according to the full PC prior, so that a sampler working on `Y`
#' only needs the standard Gaussian density as prior term. With
#' `A = sqrt(Y1^2 + Y2^2)`:
#'
#' * `B = f(0) + A^2 / (2 * lambda_v)` is the distance `f(||v||)`;
#' * `(v1, v2) = f^{-1}(B) * (Y1, Y2) / A`, replaced for `A <= eps` by its
#'   limit `(3 * pi)^(-1/4) * lambda_v^(-1/2) * (Y1, Y2)`;
#' * `kappa = -log(1 - pnorm(Y3)) / (lambda_theta * B)`;
#' * `sigma_beta = -log(1 - pnorm(Y4)) / lambda_sigma_beta`.
#'
#' `psi_reparam_vec()` is the scalar version meant for the inner loop of the
#' sampler; `psi_reparam()` wraps it and also accepts a matrix of draws.
#'
#' @param Y Numeric vector of length 4. `psi_reparam()` also accepts a matrix
#'   with 4 columns, one draw per row.
#' @param lambdas List with elements `lambda_v`, `lambda_theta` and
#'   `lambda_sigma_beta`, as returned by `pc_calibrate_all()`.
#' @param eps Threshold on `A` below which the isotropic limit is used for
#'   `(v1, v2)`.
#' @return `psi_reparam()` returns the list
#'   `list(v1, v2, kappa, sigma_beta)`, each element having one entry per
#'   draw. `psi_reparam_vec()` returns the named numeric vector
#'   `c(v1, v2, kappa, sigma_beta)`.
#' @export
psi_reparam <- function(Y, lambdas) {
  if (!is.matrix(Y)) {
    Y <- matrix(Y, nrow = 1L)
  }
  if (ncol(Y) != 4L) {
    stop("'Y' must have length 4, or be a matrix with 4 columns.")
  }
  theta <- unname(vapply(
    seq_len(nrow(Y)),
    function(i) psi_reparam_vec(Y[i, ], lambdas),
    numeric(4L)
  ))
  list(
    v1 = theta[1L, ],
    v2 = theta[2L, ],
    kappa = theta[3L, ],
    sigma_beta = theta[4L, ]
  )
}

#' @rdname psi_reparam
#' @export
psi_reparam_vec <- function(Y, lambdas, eps = 1e-8) {
  y1 <- Y[[1L]]
  y2 <- Y[[2L]]
  lambda_v <- lambdas$lambda_v

  A <- sqrt(y1 * y1 + y2 * y2)
  # f(||v||) - f(0) ~ Exp(lambda_v), since A^2 / 2 ~ Exp(1)
  excess <- A * A / (2 * lambda_v)
  B <- .pc_f0 + excess

  if (A > eps) {
    # ||v|| = f^{-1}(B) = acosh(B^2 / pi - 1/3) / 2, where
    # B^2 / pi - 1/3 = 1 + excess * (2 * f(0) + excess) / pi. Building the
    # argument from 'excess' avoids the cancellation in B^2 / pi - 4/3 and
    # keeps it >= 1 by construction.
    delta <- excess * (2 * .pc_f0 + excess) / pi
    v_scale <- 0.5 * .acosh1p(delta) / A
  } else {
    # Removable singularity: ||v|| / A -> (3 * pi)^(-1/4) * lambda_v^(-1/2)
    v_scale <- 1 / ((3 * pi)^0.25 * sqrt(lambda_v))
  }

  # -log(1 - pnorm(y)) evaluated on the log scale of the upper tail
  e_kappa <- -pnorm(Y[[3L]], lower.tail = FALSE, log.p = TRUE)
  e_sigma <- -pnorm(Y[[4L]], lower.tail = FALSE, log.p = TRUE)

  c(
    v1 = v_scale * y1,
    v2 = v_scale * y2,
    kappa = e_kappa / (lambdas$lambda_theta * B),
    sigma_beta = e_sigma / lambdas$lambda_sigma_beta
  )
}
