#' Random-effects design using restricted mean survival time
#'
#' @description Calculate sample size and shrinkage-based Method I consistency
#' probability for a treatment-minus-control restricted mean survival time
#' (RMST) difference, allowing non-proportional hazards.
#'
#' @inheritParams ssCPnorm
#' @param lambda0a,lambda0b Positive control hazard rates before and after
#'   `breakpoint`. Used only when `PW = TRUE`.
#' @param lambda1a,lambda1b Positive treatment hazard rates before and after
#'   `breakpoint`. Used only when `PW = TRUE`.
#' @param breakpoint Nonnegative, common hazard change time. Used only when `PW = TRUE`.
#' @param shape0,shape1 Positive Weibull shape parameters for control and
#'   treatment. Used only when `PW = FALSE`.
#' @param scale0,scale1 Positive Weibull scale parameters for control and
#'   treatment. Used only when `PW = FALSE`.
#' @param PW Use piecewise exponential survival (`TRUE`) or Weibull survival (`FALSE`).
#' @param unimin,unimax Uniform censoring limits, with `0 <= unimin < unimax`.
#' @param eta Positive common RMST horizon; must be below every `unimax`.
#'
#' @details
#' Supply only the survival parameters needed by `PW`. All regional model
#' parameters can be scalars or vectors of length `length(f)`.
#' Use the same time units for `eta`, `breakpoint`, censoring limits and scales.
#'
#' With `PW = TRUE`, hazards change once at `breakpoint`.
#' With `PW = FALSE`, Weibull survival is
#' \eqn{S(t)=\exp\{-(t/\mathrm{scale})^{\mathrm{shape}}\}}.
#' Scale equals the reciprocal exponential rate only when shape is one.
#' Event and censoring times are independent. Before `unimin`, censoring
#' survival equals one.
#'
#' Regional effects are RMST differences; their mean and sample standard
#' deviation define `delta` and `tau`. The mean must be positive.
#' RMST variance is integrated numerically, with
#' \eqn{\Omega_r=\sigma_{1r}^2/\ell+\sigma_{0r}^2}.
#'
#' @section Design assumptions:
#' Regional inputs may be scalars or vectors of length `length(f)`.
#' The calculation uses a normal random-effects approximation and the
#' sample-size equation in [ssCPnorm()]. Invalid inputs and designs with
#' no finite solution raise errors. CP uses the continuous design;
#' returned group counts are rounded recommendations.
#'
#' @returns A list containing `n0`, `n1`, `n` and `CP` as in [ssCPnorm()], plus:
#' \describe{
#'   \item{Dr}{Treatment-minus-control RMST difference in each region.}
#'   \item{sigma2.1r, sigma2.0r}{Treatment and control RMST variance coefficients.
#'     Divide each by its regional arm count to obtain the estimator variance.}
#' }
#' @seealso [ssCPsurvPH()], [ssCPnorm()].
#' @references Ren and Xu (2026). Consistency Assessment and Regional Sample
#' Size Calculation for MRCT Under Random Effects Model.
#' \doi{10.1002/pst.70090}.
#' @examples
#' # Piecewise exponential survival
#' ssCPsurvNPH(alpha = 0.025, beta = 0.2,
#'             lambda0a = 0.07, lambda0b = 0.04,
#'             lambda1a = c(0.02, 0.03, 0.04), lambda1b = 0.04,
#'             breakpoint = 10, unimin = 0, unimax = 240,
#'             target = 1, f = rep(1/3, 3), eta = 80)
#'
#' # Weibull survival: shape = 1 gives exponential survival
#' ssCPsurvNPH(alpha = 0.025, beta = 0.2, PW = FALSE,
#'             shape0 = 1, scale0 = 10, shape1 = 1, scale1 = 15,
#'             unimin = 20, unimax = 30, eta = 10,
#'             target = 1, f = rep(1/3, 3))
#' @export
ssCPsurvNPH<-function(alpha, beta, lambda0a, lambda0b, lambda1a, lambda1b, breakpoint,
                      shape0,scale0,shape1,scale1, PW=TRUE,
                      unimin, unimax, pi=0.5, randRatio=1, target, f, eta){

  f <- .fractions(f)
  if (length(f) < 2) stop("At least two regions are needed to estimate heterogeneity")
  .check_number(eta, "eta", 0, lower_open=TRUE)
  if (!is.logical(PW) || length(PW) != 1 || is.na(PW)) stop("PW must be TRUE or FALSE")
  unimin <- .regional(unimin, length(f), "unimin", open=FALSE)
  unimax <- .regional(unimax, length(f), "unimax")
  if (any(unimin >= unimax) || any(eta >= unimax))
    stop("Require unimin < unimax and eta < every unimax")
  if (PW) {
    .check_number(breakpoint, "breakpoint", 0)
    lambda0a <- .regional(lambda0a, length(f), "lambda0a")
    lambda0b <- .regional(lambda0b, length(f), "lambda0b")
    lambda1a <- .regional(lambda1a, length(f), "lambda1a")
    lambda1b <- .regional(lambda1b, length(f), "lambda1b")
  } else {
    shape0 <- .regional(shape0, length(f), "shape0")
    shape1 <- .regional(shape1, length(f), "shape1")
    scale0 <- .regional(scale0, length(f), "scale0")
    scale1 <- .regional(scale1, length(f), "scale1")
  }
  # Split integration at changes in the hazard or censoring model.
  split_integral <- function(fun, lower, upper, knots=numeric()) {
    points <- sort(unique(c(lower, knots[knots > lower & knots < upper], upper)))
    sum(vapply(seq_len(length(points)-1), function(i)
      integrate(fun, points[i], points[i+1], rel.tol=1e-7)$value, numeric(1)))
  }
  moments <- function(cumhaz, hazard, a, b, knots) {
    survival <- function(t) exp(-cumhaz(t))
    mean <- split_integral(survival, 0, eta, knots)
    integrand <- function(t) vapply(t, function(ti) {
      # Conditional survival avoids both inactive-branch overflow and 0/0 in tails.
      st <- survival(ti)
      if (st == 0) return(0)
      tail <- split_integral(function(u) exp(-(cumhaz(u)-cumhaz(ti))), ti, eta, knots)
      censor <- min(1, (b-ti)/(b-a))
      st * tail^2 * hazard(ti)/censor
    }, numeric(1))
    variance <- split_integral(integrand, 0, eta, c(knots, a))
    c(mean=mean, variance=variance)
  }
  arm_moments <- function(region, arm) {
    if (PW) {
      before <- if (arm == 0) lambda0a[region] else lambda1a[region]
      after <- if (arm == 0) lambda0b[region] else lambda1b[region]
      cumhaz <- function(t) before*pmin(t, breakpoint) + after*pmax(0, t-breakpoint)
      hazard <- function(t) ifelse(t <= breakpoint, before, after)
      knots <- breakpoint
    } else {
      shape <- if (arm == 0) shape0[region] else shape1[region]
      scale <- if (arm == 0) scale0[region] else scale1[region]
      cumhaz <- function(t) (t/scale)^shape
      hazard <- function(t) shape/scale*(t/scale)^(shape-1)
      knots <- numeric()
    }
    moments(cumhaz, hazard, unimin[region], unimax[region], knots)
  }
  control <- vapply(seq_along(f), arm_moments, numeric(2), arm=0)
  treatment <- vapply(seq_along(f), arm_moments, numeric(2), arm=1)
  effects <- treatment[1,]-control[1,]
  Omega <- treatment[2,]/randRatio+control[2,]
  result <- .random_design(alpha, beta, mean(effects), sd(effects), Omega,
                           pi, randRatio, target, f)
  c(list(Dr=effects, sigma2.1r=treatment[2,], sigma2.0r=control[2,]), result)
}
