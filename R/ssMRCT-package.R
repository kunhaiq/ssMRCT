#' Regional consistency and sample size for MRCTs
#'
#' @description Tools for multi-regional clinical trial (MRCT) design and
#' regional treatment-effect estimation. CP denotes a consistency probability
#' conditional on overall trial significance, not overall power or interval coverage.
#'
#' @section Choosing a function:
#' \describe{
#'   \item{\link{conProb}, \link{conProb2}}{Fixed-effects Method I CP for one or two trials.}
#'   \item{\link{regFrac}, \link{regFrac2}}{Regional fraction for a target Method I CP.
#'     The two-trial solver requires a symmetric design.}
#'   \item{\link{conProbII}, \link{conProbII2}}{Probability that all regional effects
#'     are nonnegative, for one or two trials.}
#'   \item{\link{ssCPnorm}, \link{ssCPbinary}}{Random-effects design for normal or binary endpoints.}
#'   \item{\link{ssCPsurvPH}, \link{ssCPsurvNPH}}{Random-effects design for survival endpoints.}
#'   \item{\link{CPlbd}}{Lower bound on random-effects Method I CP.}
#'   \item{\link{regATE_dr}}{Covariate-adjusted regional ATE from participant-level data.}
#' }
#'
#' @section Conventions:
#' `alpha` is one-sided. Fixed-effects functions take `power`; random-effects
#' functions take `beta = 1 - power`. Allocation ratios are treatment/control.
#' Higher treatment effects indicate benefit; PH effects are negative log hazard ratios.
#'
#' Analytical CP uses the continuous design. Returned integer sample sizes
#' are recommendations, not a re-evaluation of CP after regional rounding.
#' Binary Method II simulation instead uses the returned actual regional arm counts.
#'
#' @references
#' Qing, Kunhai, Xinru Ren, Shuping Jiang, Ping Yang, Menggang Yu, and Jin Xu. 2025.
#' \dQuote{Regional Consistency Evaluation and Sample Size Calculation under Two MRCTs.}
#' \emph{Statistics in Medicine} 44 (\enc{25–27}{25-27}): e70306.
#' \url{https://doi.org/10.1002/sim.70306}.
#'
#' Ren, Xinru, and Jin Xu. 2026.
#' \dQuote{Consistency Assessment and Regional Sample Size Calculation for MRCT
#' under Random Effects Model.} \emph{Pharmaceutical Statistics} 25 (3): e70090.
#' \url{https://doi.org/10.1002/pst.70090}.
#' @name ssMRCT-package
#' @aliases ssMRCT
#' @importFrom stats dnorm pnorm qnorm sd integrate uniroot rbinom rnorm
"_PACKAGE"
