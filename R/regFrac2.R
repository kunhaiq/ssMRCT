#' Regional fraction for two symmetric MRCTs: Method I
#'
#' @description Solve for a common regional fraction attaining the target
#' Method I consistency probability in two independent, symmetric trials.
#'
#' @inheritParams conProb
#' @param power Common overall power for both trials, with `alpha < power < 1`.
#' @param CP Target consistency probability, in `(0.5, 1]`.
#'
#' @details
#' Both trials must have the same treatment effect, arm-specific variances,
#' randomization ratio, power and regional fraction. The function solves for
#' that common fraction. It does not optimize two unequal regional fractions.
#' A target of one returns the full-region boundary.
#'
#' @returns A numeric scalar: the common regional fraction for both trials.
#' @seealso [conProb2()] for a general two-trial CP calculation.
#' @references
#' Qing, Kunhai, Xinru Ren, Shuping Jiang, Ping Yang, Menggang Yu, and Jin Xu. 2025.
#' \dQuote{Regional Consistency Evaluation and Sample Size Calculation under Two MRCTs.}
#' \emph{Statistics in Medicine} 44 (\enc{25–27}{25-27}): e70306.
#' \url{https://doi.org/10.1002/sim.70306}.
#' @examples
#' fraction <- regFrac2(alpha = 0.05, power = 0.8)
#' conProb2(alpha = 0.05, power1 = 0.8, rF1 = fraction,
#'          d1 = 1, sigmaTrt1 = 4)
#' @export
regFrac2 <- function(alpha, power, pi=0.5, CP=0.8) {
  .check_design(alpha, power, pi)
  .check_number(CP, "CP", 0.5, 1, TRUE)
  if (CP == 1) return(1)
  uniroot(function(f) .method1_cp(alpha, power, pi, f, 2L)-CP,
          c(0,1), tol=1e-8)$root
}
