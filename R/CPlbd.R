#' Random-effects consistency bound: Method I
#'
#' @description Calculate the Method I CP lower bound from Ren and Xu (2026)
#' and the regional normal-endpoint counts at its equality condition.
#'
#' @inheritParams ssCPnorm
#' @param tau_delta_r Positive ratio `tau / delta`.
#' @param fr Fractions of the regions of interest, each in `(0, 1]`.
#'   Retained for compatibility and output sizing; counts are regional and
#'   are not divided by these fractions. The entries need not sum to one.
#' @param sigma0r,sigma1r Positive control and treatment outcome standard
#'   deviations; scalars or vectors of length `length(fr)`.
#'
#' @details
#' Let \eqn{z=z_{1-\alpha}+z_{1-\beta}} and
#' \eqn{a=(\tau/\delta)^2z^2/2}.
#' If \eqn{a<1}, the function uses the bound in equation (13), with
#' regional control count at equality
#' \deqn{n_{0r}=\frac{\Omega_r}{(1/a-1)\tau^2}.}
#' If \eqn{a\ge1}, it uses equation (14), and equality counts are `NA`.
#' These counts describe the regional equality condition, not a complete
#' multi-region allocation satisfying the overall design equation.
#'
#' @returns A list with:
#' \describe{
#'   \item{CPlbd}{Scalar CP lower bound, without rounding.}
#'   \item{n0r, n1r}{Rounded regional control and treatment counts.}
#'   \item{nr}{Regional total, `n0r + n1r`.}
#' }
#' Sample-size entries are numeric `NA` when the equality condition fails.
#' @seealso [ssCPnorm()].
#' @references
#' Ren, Xinru, and Jin Xu. 2026.
#' \dQuote{Consistency Assessment and Regional Sample Size Calculation for MRCT
#' under Random Effects Model.} \emph{Pharmaceutical Statistics} 25 (3): e70090.
#' \doi{10.1002/pst.70090}.
#' @examples
#' CPlbd(alpha = 0.05, beta = 0.2, tau_delta_r = 0.4, delta = 0.25,
#'       sigma0r = 1, sigma1r = 1, fr = c(0.1, 0.5))
#' @export
CPlbd <- function(alpha, beta, tau_delta_r, delta, sigma0r, sigma1r, fr,
                  pi = 0.5, randRatio = 1) {
  .check_number(beta, "beta", 0, 1, TRUE, TRUE)
  .check_design(alpha, 1 - beta, pi)
  .check_number(tau_delta_r, "tau_delta_r", 0, lower_open = TRUE)
  .check_number(delta, "delta", 0, lower_open = TRUE)
  .check_number(randRatio, "randRatio", 0, lower_open = TRUE)
  .check_number(fr, "fr", 0, 1, TRUE, FALSE, FALSE)
  sigma0r <- .regional(sigma0r, length(fr), "sigma0r")
  sigma1r <- .regional(sigma1r, length(fr), "sigma1r")
  z <- qnorm(1 - alpha) + qnorm(1 - beta)
  h <- tau_delta_r * z

  if (h < sqrt(2)) {
    # Equation (13); this form avoids cancellation when heterogeneity is small.
    slope <- 2 * (1 - pi) / h^2
    Omega <- sigma0r^2 + sigma1r^2 / randRatio
    n0r <- ceiling(Omega * (z / delta)^2 / (2 - h^2))
    n1r <- ceiling(n0r * randRatio)
    nr <- n0r + n1r
  } else {
    # Equation (14) has no finite regional equality counts.
    slope <- (1 - pi) / sqrt(h^2 - 1)
    n0r <- n1r <- nr <- rep(NA_real_, length(fr))
  }
  integrand <- function(x) pnorm(slope * (x + z)) * dnorm(x)
  cpmin <- integrate(integrand, lower = qnorm(beta), upper = Inf)$value / (1 - beta)
  list(CPlbd = cpmin, n0r = n0r, n1r = n1r, nr = nr)
}
