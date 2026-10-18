#' Random-effects design for a binary endpoint: Method I
#'
#' @description Calculate overall sample size and shrinkage-based Method I
#' consistency probability from regional response probabilities.
#'
#' @inheritParams ssCPnorm
#' @param p0r,p1r Control and treatment response probabilities, strictly in
#'   `(0, 1)`, as scalars or vectors with one value per region.
#'
#' @details
#' Regional effects are `p1r - p0r`. Their unweighted mean and sample
#' standard deviation define `delta` and `tau`; the mean must be positive.
#' The sampling variance component is
#' \deqn{\Omega_r = p_{1r}(1-p_{1r})/\ell+p_{0r}(1-p_{0r}).}
#' This function uses a normal random-effects approximation, not binomial
#' simulation. For fixed-effects Method II simulation, use [conProbII()].
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
#' @seealso [ssCPnorm()], [CPlbd()].
#' @references
#' Ren, Xinru, and Jin Xu. 2026.
#' \dQuote{Consistency Assessment and Regional Sample Size Calculation for MRCT
#' under Random Effects Model.} \emph{Pharmaceutical Statistics} 25 (3): e70090.
#' \url{https://doi.org/10.1002/pst.70090}.
#' @examples
#' ssCPbinary(alpha = 0.025, beta = 0.2, p0r = 0.3,
#'            p1r = c(0.6, 0.55, 0.5), target = 1, f = rep(1/3, 3))
#' @export
ssCPbinary <- function(alpha, beta, p0r, p1r, pi=0.5, randRatio=1, target, f) {
  f <- .fractions(f)
  if (length(f) < 2) stop("At least two regions are needed to estimate heterogeneity")
  p0r <- .regional(p0r, length(f), "p0r", upper=1)
  p1r <- .regional(p1r, length(f), "p1r", upper=1)
  if (any(p0r == 1 | p1r == 1)) stop("Response probabilities must be strictly below 1")
  effects <- p1r-p0r
  .random_design(alpha, beta, mean(effects), sd(effects),
    p1r*(1-p1r)/randRatio+p0r*(1-p0r), pi, randRatio, target, f)
}
