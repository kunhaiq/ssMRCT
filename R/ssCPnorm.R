#' Random-effects design for a normal endpoint: Method I
#'
#' @description Calculate overall sample size and conditional Method I
#' consistency probability for an empirical shrinkage regional estimator.
#'
#' @param alpha One-sided significance level, between 0 and 0.5.
#' @param beta Type II error probability; `alpha < 1 - beta < 1`.
#' @param delta Positive mean of the regional treatment-effect distribution.
#' @param tau Nonnegative standard deviation of regional treatment effects.
#' @param sigma1r,sigma0r Positive treatment and control outcome standard
#'   deviations, as scalars or vectors with one value per region.
#' @param pi Fraction of the overall effect to preserve; `0 <= pi < 1`.
#' @param randRatio Treatment-to-control allocation ratio, strictly positive.
#' @param target Index of the region of interest in `f`.
#' @param f Positive regional fractions summing to one. At least two regions are required.
#'
#' @details
#' Let \eqn{\ell} be `randRatio`, \eqn{z = z_{1-\alpha}+z_{1-\beta}}
#' and \eqn{\Omega_r = \sigma_{1r}^2/\ell+\sigma_{0r}^2}.
#' The continuous control count solves
#' \deqn{\sum_r \left(\tau^2+\frac{\Omega_r}{n_0 f_r}\right)^{-1}
#'   = \frac{z^2}{\delta^2}.}
#' A finite solution requires \eqn{\tau z/\delta < \sqrt{R}},
#' where \eqn{R} is the number of regions.
#'
#' The target-region criterion is
#' \deqn{\Pr(\widetilde D_r \ge \pi\widetilde D
#'   \mid T > z_{1-\alpha}).}
#' It uses the shrinkage estimator, not the raw regional mean difference.
#' At `tau = 0`, regional and overall shrinkage estimates coincide, so CP is one.
#'
#' @section Design assumptions:
#' Regional inputs may be scalars or vectors of length `length(f)`.
#' The calculation uses a normal random-effects approximation and the
#' sample-size equation above. Invalid inputs and designs with
#' no finite solution raise errors. CP uses the continuous design;
#' returned group counts are rounded recommendations.
#'
#' @returns A list with:
#' \describe{
#'   \item{n0, n1}{Recommended overall control and treatment counts.}
#'   \item{n}{Total count, `n0 + n1`.}
#'   \item{CP}{Conditional consistency probability, without rounding.}
#' }
#' @seealso [ssCPbinary()], [ssCPsurvPH()], [ssCPsurvNPH()], [CPlbd()].
#' @references
#' Ren, Xinru, and Jin Xu. 2026.
#' \dQuote{Consistency Assessment and Regional Sample Size Calculation for MRCT
#' under Random Effects Model.} \emph{Pharmaceutical Statistics} 25 (3): e70090.
#' \url{https://doi.org/10.1002/pst.70090}.
#' @examples
#' ssCPnorm(alpha = 0.025, beta = 0.2, delta = 0.25, tau = 0.1,
#'          sigma1r = 1, sigma0r = 1, target = 1, f = rep(1/3, 3))
#' @export
ssCPnorm <- function(alpha, beta, delta, tau, sigma1r, sigma0r,
                     pi=0.5, randRatio=1, target, f) {
  f <- .fractions(f)
  sigma1r <- .regional(sigma1r, length(f), "sigma1r")
  sigma0r <- .regional(sigma0r, length(f), "sigma0r")
  .random_design(alpha, beta, delta, tau, sigma1r^2/randRatio+sigma0r^2,
                 pi, randRatio, target, f)
}
