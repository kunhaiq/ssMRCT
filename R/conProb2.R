#' Pooled consistency probability for two MRCTs: Method I
#'
#' @description Calculate Method I consistency probability from pooled data,
#' conditional on both independent trials being significant.
#'
#' @param alpha One-sided significance level, between 0 and 0.5.
#' @param power1,power2 Overall power in studies 1 and 2; each lies in `(alpha, 1)`.
#' @param d1,d2 Positive treatment-minus-control mean differences.
#' @param sigmaTrt1,sigmaTrt2 Positive treatment-arm outcome standard deviations.
#' @param sigmaCtrl1,sigmaCtrl2 Positive control-arm outcome standard deviations.
#' @param randRatio1,randRatio2 Positive treatment-to-control allocation ratios.
#' @param muCtrl1,muCtrl2 Control-arm means; needed when pooling weights differ.
#' @param pi Fraction of the overall effect to preserve; `0 <= pi < 1`.
#' @param rF1,rF2 Fractions in the target region in studies 1 and 2, in `(0, 1]`.
#'
#' @details
#' Treatment and control outcomes are pooled separately using their sample
#' sizes. Overall and regional pooling weights can differ.
#' \deqn{\Pr(D_{r,\mathrm{pool}} \ge \pi D_{\mathrm{pool}}
#'   \mid T_1 > z_{1-\alpha}, T_2 > z_{1-\alpha}).}
#' The conditional normal probability is integrated over both significant
#' study results. Study 2 defaults to the corresponding study 1 inputs.
#'
#' @section Sample sizes:
#' CP uses the continuous design at the requested power. The recommended
#' control count is rounded up first, followed by the treatment count.
#' Returned counts do not re-evaluate CP after integer regional allocation.
#'
#' @returns A list with `CP`, `N1` and `N2`: the conditional consistency
#' probability and recommended total sample sizes in studies 1 and 2.
#' @seealso [conProb()], [regFrac2()].
#' @references Qing et al. (2025). Regional consistency evaluation and sample
#' size calculation under two MRCTs. \doi{10.1002/sim.70306}.
#' The implementation uses the corrected pooling weights and covariance structure.
#' @examples
#' conProb2(alpha = 0.05, power1 = 0.8, power2 = 0.9,
#'          rF1 = 0.141, d1 = 1, sigmaTrt1 = 4)
#' @export
conProb2 <- function(alpha,
                     power1,
                     power2 = power1,
                     pi = 0.5,
                     rF1,
                     rF2 = rF1,
                     d1,
                     d2 = d1,
                     muCtrl1 = 0,
                     muCtrl2 = muCtrl1,
                     sigmaTrt1,
                     sigmaCtrl1 = sigmaTrt1,
                     sigmaTrt2 = sigmaTrt1,
                     sigmaCtrl2 = sigmaTrt2,
                     randRatio1 = 1,
                     randRatio2 = randRatio1){

  .check_design(alpha, power1, pi)
  .check_design(alpha, power2, pi)
  .check_number(rF1, "rF1", 0, 1, TRUE)
  .check_number(rF2, "rF2", 0, 1, TRUE)
  .check_number(muCtrl1, "muCtrl1")
  .check_number(muCtrl2, "muCtrl2")
  des1 <- .fixed_design(alpha, power1, d1, sigmaTrt1, sigmaCtrl1, randRatio1)
  des2 <- .fixed_design(alpha, power2, d2, sigmaTrt2, sigmaCtrl2, randRatio2)
  N1 <- des1$nTrt + des1$nCtrl
  N2 <- des2$nTrt + des2$nCtrl
  if (rF1 == 1 && rF2 == 1 && randRatio1 == randRatio2) {
    return(list(CP=1, N1=N1, N2=N2))
  }

  beta1 <- 1 - power1
  beta2 <- 1 - power2
  mu_c1 <- muCtrl1
  mu_c2 <- muCtrl2
  r1 <- randRatio1
  r2 <- randRatio2
  fk1 <- rF1
  fk2 <- rF2
  sigma_t1 <- sigmaTrt1; sigma_c1 <- sigmaCtrl1
  sigma_t2 <- sigmaTrt2; sigma_c2 <- sigmaCtrl2

  N_t1 <- des1$nt
  N_c1 <- des1$nc
  N_t2 <- des2$nt
  N_c2 <- des2$nc

  Nk_t1 <- fk1 * N_t1; Nk_c1 <- fk1 * N_c1
  Nk_t2 <- fk2 * N_t2; Nk_c2 <- fk2 * N_c2

  mu_t1 <- mu_c1 + d1
  mu_t2 <- mu_c2 + d2

  Denom_T <- N_t1 + N_t2
  Denom_C <- N_c1 + N_c2

  W_t1 <- N_t1 / Denom_T; W_t2 <- N_t2 / Denom_T
  W_c1 <- N_c1 / Denom_C; W_c2 <- N_c2 / Denom_C

  # Regional and overall pooling weights can differ.
  Denom_k_T <- Nk_t1 + Nk_t2
  Denom_k_C <- Nk_c1 + Nk_c2

  wk_t1 <- Nk_t1 / Denom_k_T; wk_t2 <- Nk_t2 / Denom_k_T
  wk_c1 <- Nk_c1 / Denom_k_C; wk_c2 <- Nk_c2 / Denom_k_C

  vY_t1 <- sigma_t1^2 / N_t1; vY_c1 <- sigma_c1^2 / N_c1
  vY_t2 <- sigma_t2^2 / N_t2; vY_c2 <- sigma_c2^2 / N_c2

  sig2_d1 <- vY_t1 + vY_c1
  sig2_d2 <- vY_t2 + vY_c2
  sig_d1 <- sqrt(sig2_d1)
  sig_d2 <- sqrt(sig2_d2)

  lam_1 <- wk_t1 * vY_t1 + wk_c1 * vY_c1
  lam_2 <- wk_t2 * vY_t2 + wk_c2 * vY_c2

  Lam_1 <- W_t1 * vY_t1 + W_c1 * vY_c1
  Lam_2 <- W_t2 * vY_t2 + W_c2 * vY_c2

  C1 <- lam_1 - pi * Lam_1
  C2 <- lam_2 - pi * Lam_2

  Term_Region <- (wk_t1 * mu_t1 - wk_c1 * mu_c1) + (wk_t2 * mu_t2 - wk_c2 * mu_c2)
  Term_Global <- (W_t1 * mu_t1 - W_c1 * mu_c1) + (W_t2 * mu_t2 - W_c2 * mu_c2)

  EL <- Term_Region - pi * Term_Global

  calc_var_term <- function(wk, W, vY, fk) {
    (wk^2 * vY / fk) + (pi^2 * W^2 * vY) - (2 * pi * wk * W * vY)
  }

  varL <- 0
  varL <- varL + calc_var_term(wk_t1, W_t1, vY_t1, fk1)
  varL <- varL + calc_var_term(wk_c1, W_c1, vY_c1, fk1)
  varL <- varL + calc_var_term(wk_t2, W_t2, vY_t2, fk2)
  varL <- varL + calc_var_term(wk_c2, W_c2, vY_c2, fk2)

  # Conditional variance of L = regional effect - pi * overall effect.
  sig2_LD <- varL - (C1^2 / sig2_d1) - (C2^2 / sig2_d2)

  if (sig2_LD < -1e-10 * varL) stop("Invalid conditional variance")
  if (abs(sig2_LD) <= 1e-10 * varL) sig2_LD <- 0
  sig2_LD <- max(0, sig2_LD)
  sig_LD <- sqrt(sig2_LD)

  z_b1 <- qnorm(1 - beta1)
  z_b2 <- qnorm(1 - beta2)

  # Average the conditional probability over two significant study results.
  integrand_func <- function(v) {
    sapply(v, function(vi) {
      inner_func <- function(u) {
        mu_val <- EL + (C1 / sig_d1) * u + (C2 / sig_d2) * vi
        .normal_nonnegative(mu_val, sig_LD) * dnorm(u)
      }
      res_u <- integrate(inner_func, lower = -z_b1, upper = Inf)$value
      return(res_u * dnorm(vi))
    })
  }

  res_total <- integrate(integrand_func, lower = -z_b2, upper = Inf)$value

  final_cp <- res_total / ((1 - beta1) * (1 - beta2))

  return(list(CP = as.numeric(final_cp), N1 = N1, N2 = N2))
}
