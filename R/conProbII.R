#' Consistency probability for all regions in one MRCT
#'
#' @description Calculate Method II consistency probability: all regional
#' effect estimates are nonnegative, conditional on overall significance.
#'
#' @inheritParams conProb
#' @param rF Positive fractions for all regions, summing to one.
#' @param d,sigmaTrt,sigmaCtrl Continuous-endpoint effect and standard deviations;
#'   see [conProb()]. Ignored for binary endpoints.
#' @param pTrt,pCtrl Treatment and control response probabilities for a binary
#'   endpoint, strictly between zero and one, with `pTrt > pCtrl`.
#' @param responseType Endpoint type: `"continuous"` or `"binary"`.
#' @param B Positive integer number of replicates; used only for binary endpoints.
#'
#' @details
#' \deqn{\mathrm{CP} = \Pr(D_r \ge 0 \text{ for all }r
#'   \mid T > z_{1-\alpha}).}
#' Continuous endpoints use a joint normal probability, retaining dependence
#' between regional and overall effects. Numerical integration can vary
#' slightly with the random seed.
#'
#' @section Sample sizes:
#' CP uses the continuous design at the requested power. The recommended
#' control count is rounded up first, followed by the treatment count.
#' Returned counts do not re-evaluate CP after integer regional allocation.
#'
#' @section Binary calculation:
#' The function simulates regional binomial counts with separately rounded
#' arm sizes. CP is the number of consistent overall successes divided by
#' the number of overall successes. Zero regional effects count as consistent.
#' Wald tests use the actual arm counts. With zero estimated variance, only
#' a positive estimated effect rejects the null.
#'
#' No overall successes gives `NA` for CP and its Monte Carlo standard error,
#' with a warning. Set the random seed for reproducibility.
#'
#' @returns A list with `CP` and total sample size `N`.
#' For binary endpoints, `N` is the actual simulated total and the list also contains:
#' \describe{
#'   \item{nTrt, nCtrl}{Actual overall treatment and control counts.}
#'   \item{nTrtReg, nCtrlReg}{Actual regional treatment and control counts.}
#'   \item{achievedPower}{Simulated probability of the conditioning event.}
#'   \item{MCSE}{Monte Carlo standard error of CP.}
#'   \item{successes}{Number of replicates meeting the conditioning event.}
#'   \item{consistentSuccesses}{Number also meeting regional consistency.}
#'   \item{B}{Number of simulation replicates.}
#' }
#' @seealso [conProbII2()], [conProb()].
#' @references Qing et al. (2025). Regional consistency evaluation and sample
#' size calculation under two MRCTs. \doi{10.1002/sim.70306}.
#' The implementation uses the corrected pooling weights and covariance structure.
#' @examples
#' set.seed(42)
#' conProbII(alpha = 0.05, power = 0.8, rF = rep(1/3, 3),
#'           d = 1, sigmaTrt = 4)
#' conProbII(alpha = 0.05, power = 0.8, rF = rep(1/3, 3),
#'           pTrt = 0.8, pCtrl = 0.7, responseType = "binary", B = 10000)
#' @importFrom mvtnorm pmvnorm GenzBretz
#' @export
conProbII <- function(alpha, power, rF, d, sigmaTrt,
                      sigmaCtrl=sigmaTrt, pTrt, pCtrl, randRatio=1,
                      responseType=c("continuous", "binary"), B=100000) {
  responseType <- match.arg(responseType)
  .check_design(alpha, power)
  rF <- .fractions(rF, "rF")
  if (responseType == "binary") {
    des <- .binary_design(alpha, power, pTrt, pCtrl, randRatio, rF)
    ans <- .binary_cp(list(des), alpha, B)
    return(c(ans, list(N=sum(des$nTrtReg)+sum(des$nCtrlReg),
                       nTrt=sum(des$nTrtReg), nCtrl=sum(des$nCtrlReg),
                       nTrtReg=des$nTrtReg, nCtrlReg=des$nCtrlReg)))
  }
  des <- .fixed_design(alpha, power, d, sigmaTrt, sigmaCtrl, randRatio)
  k <- length(rF)
  if (k == 1L) return(list(CP=1, N=des$nTrt+des$nCtrl))
  # Integrate the joint event; cross-covariances preserve the dependence on D.
  Sigma <- diag(c(1/rF, 1))
  Sigma[seq_len(k), k+1] <- Sigma[k+1, seq_len(k)] <- 1
  z <- qnorm(1-alpha)+qnorm(power)
  prob <- mvtnorm::pmvnorm(lower=c(rep(0,k), -qnorm(power)),
    upper=rep(Inf,k+1), mean=c(rep(z,k),0), sigma=Sigma,
    algorithm=mvtnorm::GenzBretz(maxpts=1000000, abseps=1e-6))
  list(CP=as.numeric(prob)/power, N=des$nTrt+des$nCtrl)
}
