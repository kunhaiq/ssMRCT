#' Regional fraction for one MRCT: Method I
#'
#' @description Find the smallest continuous regional fraction attaining a
#' target consistency probability under the fixed-effects model.
#'
#' @inheritParams conProb
#' @param CP Target consistency probability, in `(0.5, 1]`.
#'
#' @details
#' Solves `conProb(...)$CP = CP` for `rF`. A target of one returns `rF = 1`.
#' Fractions are not converted to integer regional arm counts.
#'
#' @section Sample sizes:
#' CP uses the continuous design at the requested power. The recommended
#' control count is rounded up first, followed by the treatment count.
#' Returned counts do not re-evaluate CP after integer regional allocation.
#'
#' @returns A list with `rF` (regional fraction) and `N` (recommended total count).
#' @seealso [conProb()], [regFrac2()].
#' @examples
#' regFrac(alpha = 0.05, power = 0.8, CP = 0.8)
#' @export
regFrac <- function(alpha, power, pi=0.5, CP=0.8, d=1,
                    sigmaTrt=4, sigmaCtrl=sigmaTrt, randRatio=1) {
  .check_design(alpha, power, pi)
  .check_number(CP, "CP", 0.5, 1, TRUE)
  des <- .fixed_design(alpha, power, d, sigmaTrt, sigmaCtrl, randRatio)
  rF <- if (CP == 1) 1 else uniroot(function(f)
    .method1_cp(alpha, power, pi, f)-CP, c(0,1), tol=1e-8)$root
  list(rF=rF, N=des$nTrt+des$nCtrl)
}
