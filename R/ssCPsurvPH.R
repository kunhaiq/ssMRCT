#' Random-effects design for a PH endpoint: Method I
#'
#' @description Calculate sample size and shrinkage-based Method I consistency
#' probability for exponential survival with a common fixed follow-up time.
#'
#' @inheritParams ssCPnorm
#' @param lambda0r Positive control hazard rates, as a scalar or regional vector.
#' @param hr Positive treatment-to-control hazard ratios, as a scalar or
#'   regional vector. The mean of `-log(hr)` must be positive.
#' @param L Positive follow-up duration, in the same time units as the hazard rates.
#'
#' @details
#' All participants are followed for `L` time units. Administrative censoring
#' at `L` is the only censoring mechanism. Regional effects are `-log(hr)`;
#' their mean and sample standard deviation define `delta` and `tau`.
#' With event probabilities \eqn{P_{kr}=1-\exp(-\lambda_{kr}L)},
#' \deqn{\Omega_r=\frac{(\ell+1)^2}{\ell(P_{0r}+\ell P_{1r})}.}
#' For delayed effects or crossing survival curves, use [ssCPsurvNPH()].
#'
#' @section Design assumptions:
#' Regional inputs may be scalars or vectors of length `length(f)`.
#' The calculation uses a normal random-effects approximation and the
#' sample-size equation in [ssCPnorm()]. Invalid inputs and designs with
#' no finite solution raise errors. CP uses the continuous design;
#' returned group counts are rounded recommendations.
#'
#' @returns A list with:
#' \describe{
#'   \item{n0, n1}{Recommended overall control and treatment counts.}
#'   \item{n}{Total count, `n0 + n1`.}
#'   \item{CP}{Conditional consistency probability, without rounding.}
#' }
#' @seealso [ssCPnorm()], [ssCPsurvNPH()].
#' @references
#' Ren, Xinru, and Jin Xu. 2026.
#' \dQuote{Consistency Assessment and Regional Sample Size Calculation for MRCT
#' under Random Effects Model.} \emph{Pharmaceutical Statistics} 25 (3): e70090.
#' \url{https://doi.org/10.1002/pst.70090}.
#' @examples
#' ssCPsurvPH(alpha = 0.025, beta = 0.2, lambda0r = 0.05,
#'            hr = c(0.7, 0.6, 0.4), target = 1, f = rep(1/3, 3), L = 36)
#' @export
ssCPsurvPH <- function(alpha, beta, lambda0r, hr, pi=0.5, randRatio=1, target, f, L) {
  f <- .fractions(f)
  if (length(f) < 2) stop("At least two regions are needed to estimate heterogeneity")
  lambda0r <- .regional(lambda0r, length(f), "lambda0r")
  hr <- .regional(hr, length(f), "hr")
  .check_number(L, "L", 0, lower_open=TRUE)
  .check_number(randRatio, "randRatio", 0, lower_open=TRUE)
  effects <- -log(hr)
  P0 <- -expm1(-lambda0r*L)
  P1 <- -expm1(-lambda0r*hr*L)
  Omega <- (randRatio+1)^2/(randRatio*(P0+randRatio*P1))
  .random_design(alpha, beta, mean(effects), sd(effects), Omega, pi, randRatio, target, f)
}
