#' Pooled consistency probability for all regions in two MRCTs
#'
#' @description Calculate Method II consistency probability from pooled
#' regional data, conditional on both independent trials being significant.
#'
#' @inheritParams conProb2
#' @inheritParams conProbII
#' @param rF1,rF2 Positive regional fraction vectors, each summing to one
#'   and listing the same regions in the same order.
#' @param d1,d2 Positive continuous-endpoint effects; ignored for binary endpoints.
#' @param pTrt1,pTrt2 Treatment response probabilities for binary endpoints,
#'   strictly between zero and one.
#' @param pCtrl1,pCtrl2 Control response probabilities for binary endpoints,
#'   strictly between zero and the corresponding treatment probability.
#' @param muCtrl1,muCtrl2 Control means for continuous endpoints only.
#'
#' @details
#' Consistency requires every pooled regional effect to be nonnegative.
#' Treatment and control outcomes are pooled separately. The conditioning
#' event is significance in both trials, not in a single pooled test.
#' For continuous endpoints, a joint normal calculation accounts for
#' regional/overall dependence and different control means. Set the seed
#' for reproducible numerical integration or simulation.
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
#' @returns A list with `CP`, `N1` and `N2`. For binary endpoints, `N1` and
#' `N2` are the actual simulated totals and additional entries are:
#' \describe{
#'   \item{nTrtReg1, nCtrlReg1}{Regional arm counts in study 1.}
#'   \item{nTrtReg2, nCtrlReg2}{Regional arm counts in study 2.}
#'   \item{achievedPower}{Simulated probability of the conditioning event.}
#'   \item{MCSE}{Monte Carlo standard error of CP.}
#'   \item{successes}{Number of replicates meeting the conditioning event.}
#'   \item{consistentSuccesses}{Number also meeting regional consistency.}
#'   \item{B}{Number of simulation replicates.}
#' }
#' Here `achievedPower` is the simulated probability that both studies succeed.
#' @seealso [conProbII()], [conProb2()].
#' @references Qing et al. (2025). Regional consistency evaluation and sample
#' size calculation under two MRCTs. \doi{10.1002/sim.70306}.
#' The implementation uses the corrected pooling weights and covariance structure.
#' @examples
#' set.seed(42)
#' conProbII2(alpha = 0.05, power1 = 0.8, rF1 = rep(1/3, 3),
#'            d1 = 1, sigmaTrt1 = 4)
#' @export
conProbII2 <- function(alpha,
                       power1,
                       power2 = power1,
                       rF1,
                       rF2 = rF1,
                       d1,
                       d2 = d1,
                       sigmaTrt1,
                       sigmaCtrl1 = sigmaTrt1,
                       sigmaTrt2 = sigmaTrt1,
                       sigmaCtrl2 = sigmaTrt2,
                       pTrt1,
                       pCtrl1,
                       pTrt2 = pTrt1,
                       pCtrl2 = pCtrl1,
                       randRatio1 = 1,
                       randRatio2 = randRatio1,
                       responseType = c("continuous","binary"),
                       B = 100000,
                       muCtrl1 = 0, muCtrl2 = muCtrl1){

  responseType <- match.arg(responseType)
  .check_design(alpha, power1)
  .check_design(alpha, power2)
  rF1 <- .fractions(rF1, "rF1")
  rF2 <- .fractions(rF2, "rF2")
  if (length(rF1) != length(rF2)) stop("rF1 and rF2 must have equal lengths")
  if (responseType == "binary") {
    des1 <- .binary_design(alpha, power1, pTrt1, pCtrl1, randRatio1, rF1)
    des2 <- .binary_design(alpha, power2, pTrt2, pCtrl2, randRatio2, rF2)
    ans <- .binary_cp(list(des1, des2), alpha, B)
    return(c(ans, list(N1=sum(des1$nTrtReg)+sum(des1$nCtrlReg),
                       N2=sum(des2$nTrtReg)+sum(des2$nCtrlReg),
                       nTrtReg1=des1$nTrtReg, nCtrlReg1=des1$nCtrlReg,
                       nTrtReg2=des2$nTrtReg, nCtrlReg2=des2$nCtrlReg)))
  }
  .check_number(muCtrl1, "muCtrl1")
  .check_number(muCtrl2, "muCtrl2")
  des1 <- .fixed_design(alpha, power1, d1, sigmaTrt1, sigmaCtrl1, randRatio1)
  des2 <- .fixed_design(alpha, power2, d2, sigmaTrt2, sigmaCtrl2, randRatio2)
  if (responseType == "continuous"){
    N1 <- des1$nc + des1$nt
    N2 <- des2$nc + des2$nt

    Nt1 <- N1 * randRatio1 / (randRatio1 + 1)
    Nc1 <- N1 * 1 / (randRatio1 + 1)
    Nt2 <- N2 * randRatio2 / (randRatio2 + 1)
    Nc2 <- N2 * 1 / (randRatio2 + 1)

    var_Ybar_t1 <- sigmaTrt1^2 / Nt1
    var_Ybar_c1 <- sigmaCtrl1^2 / Nc1
    var_Ybar_t2 <- sigmaTrt2^2 / Nt2
    var_Ybar_c2 <- sigmaCtrl2^2 / Nc2

    sigmad1_sq <- var_Ybar_t1 + var_Ybar_c1
    sigmad2_sq <- var_Ybar_t2 + var_Ybar_c2

    K <- length(rF1)

    Nt1_k <- rF1 * Nt1
    Nc1_k <- rF1 * Nc1
    Nt2_k <- rF2 * Nt2
    Nc2_k <- rF2 * Nc2

    denom_t <- Nt1_k + Nt2_k
    denom_c <- Nc1_k + Nc2_k

    w_t1 <- Nt1_k / denom_t
    w_t2 <- Nt2_k / denom_t
    w_c1 <- Nc1_k / denom_c
    w_c2 <- Nc2_k / denom_c

    lambda1 <- w_t1 * var_Ybar_t1 + w_c1 * var_Ybar_c1
    lambda2 <- w_t2 * var_Ybar_t2 + w_c2 * var_Ybar_c2

    var_Dk_pool <- (w_t1^2 * sigmaTrt1^2 / Nt1_k) + (w_c1^2 * sigmaCtrl1^2 / Nc1_k) +
      (w_t2^2 * sigmaTrt2^2 / Nt2_k) + (w_c2^2 * sigmaCtrl2^2 / Nc2_k)
    Sigma_DD <- diag(var_Dk_pool, nrow=K, ncol=K)

    Cov_DZ1 <- lambda1 / sqrt(sigmad1_sq) # Vector length K
    Cov_DZ2 <- lambda2 / sqrt(sigmad2_sq) # Vector length K

    # Joint covariance retains dependence between regional and overall effects.
    Sigma_total <- matrix(0, nrow = K+2, ncol = K+2)
    Sigma_total[1:K, 1:K] <- Sigma_DD
    Sigma_total[(K+1), (K+1)] <- 1
    Sigma_total[(K+2), (K+2)] <- 1

    Sigma_total[1:K, K+1] <- Cov_DZ1
    Sigma_total[K+1, 1:K] <- Cov_DZ1
    Sigma_total[1:K, K+2] <- Cov_DZ2
    Sigma_total[K+2, 1:K] <- Cov_DZ2

    # Unequal arm weights require the control means as well as treatment effects.
    mean_D <- w_t1 * (muCtrl1+d1) - w_c1*muCtrl1 +
      w_t2 * (muCtrl2+d2) - w_c2*muCtrl2
    Mean_total <- c(mean_D, 0, 0)

    lower_limits <- c(rep(0, K), -qnorm(power1), -qnorm(power2))
    upper_limits <- rep(Inf, K+2)

    prob_joint <- mvtnorm::pmvnorm(lower = lower_limits,
                                   upper = upper_limits,
                                   mean = Mean_total,
                                   sigma = Sigma_total,
                                   algorithm=mvtnorm::GenzBretz(maxpts=1000000, abseps=1e-6))

    CP <- as.numeric(prob_joint) / (power1 * power2)

    return(list(CP=CP, N1=des1$nTrt+des1$nCtrl, N2=des2$nTrt+des2$nCtrl))
  }
}
