#' Consistency probability for one MRCT: Method I
#'
#' @description Calculate the probability that a regional treatment effect
#' preserves a specified fraction of the overall effect, conditional on
#' overall significance, under a fixed-effects model.
#'
#' @param alpha One-sided significance level, between 0 and 0.5.
#' @param power Target overall power, with `alpha < power < 1`.
#' @param pi Fraction of the overall effect to preserve; `0 <= pi < 1`.
#' @param rF Fraction of participants in the region of interest, in `(0, 1]`.
#' @param d Positive treatment-minus-control mean difference.
#' @param sigmaTrt Positive outcome standard deviation in the treatment arm.
#' @param sigmaCtrl Positive outcome standard deviation in the control arm.
#' @param randRatio Treatment-to-control allocation ratio, strictly positive.
#'
#' @details
#' With regional effect estimate \eqn{D_r} and overall estimate \eqn{D},
#' \deqn{\mathrm{CP} = \Pr(D_r \ge \pi D \mid T > z_{1-\alpha}).}
#' The normal approximation is integrated numerically. CP depends on
#' `alpha`, `power`, `pi` and `rF`; the other inputs determine sample size.
#'
#' @section Sample sizes:
#' CP uses the continuous design at the requested power. The recommended
#' control count is rounded up first, followed by the treatment count.
#' Returned counts do not re-evaluate CP after integer regional allocation.
#'
#' @returns A list with `CP` (conditional consistency probability) and `N`
#' (recommended total sample size).
#' @seealso [regFrac()] to solve for a regional fraction; [conProb2()] for two trials.
#' @references
#' Qing, Kunhai, Xinru Ren, Shuping Jiang, Ping Yang, Menggang Yu, and Jin Xu. 2025.
#' \dQuote{Regional Consistency Evaluation and Sample Size Calculation under Two MRCTs.}
#' \emph{Statistics in Medicine} 44 (\enc{25–27}{25-27}): e70306.
#' \doi{10.1002/sim.70306}.
#' @examples
#' conProb(alpha = 0.05, power = 0.8, rF = 0.271,
#'         d = 1, sigmaTrt = 4)
#' @export
conProb <- function(alpha, power, pi=0.5, rF, d, sigmaTrt,
                    sigmaCtrl=sigmaTrt, randRatio=1) {
  .check_design(alpha, power, pi)
  .check_number(rF, "rF", 0, 1, TRUE)
  des <- .fixed_design(alpha, power, d, sigmaTrt, sigmaCtrl, randRatio)
  list(CP=.method1_cp(alpha, power, pi, rF), N=des$nTrt+des$nCtrl)
}
