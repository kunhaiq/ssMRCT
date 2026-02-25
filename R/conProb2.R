#' @name conProb2
#'
#' @title Consistency probability for two MRCTs via the extended criterion Method I of MHLW (2007) (conditional version) under fixed effects model
#'
#' @description Calculate the consistency probability for two MRCTs via the extended criterion Method I of MHLW (2007) (conditional version) under fixed effects model.
#'
#' @param alpha The Type I error.
#' @param power1 Power for MRCT 1.
#' @param power2 Power for MRCT 2. Defaults to \code{power1}.
#' @param pi The threshold ratio in the extended criterion Method I of MHLW (2007) (conditional version). Defaults to 0.5.
#' @param rF1 The regional fraction for MRCT 1.
#' @param rF2 The regional fraction for MRCT 2. Defaults to \code{rF1}.
#' @param d1 The true mean of difference of response for MRCT 1.
#' @param d2 The true mean of difference of response for MRCT 2. Defaults to \code{d1}.
#' @param muCtrl1 The true mean of response in the control group for MRCT 1. Defaults to 0.
#' @param muCtrl2 The true mean of response in the control group for MRCT 2. Defaults to \code{muCtrl1}.
#' @param sigmaTrt1 The standard deviation of response in the treatment group for MRCT 1.
#' @param sigmaCtrl1 The standard deviation of response in the control group for MRCT 1. Defaults to \code{sigmaTrt1}.
#' @param sigmaTrt2 The standard deviation of response in the treatment group for MRCT 2. Defaults to \code{sigmaTrt1}.
#' @param sigmaCtrl2 The standard deviation of response in the control group for MRCT 2. Defaults to \code{sigmaTrt2}.
#' @param randRatio1 The randomization ratio between the treatment group and control group for MRCT 1. Defaults to 1.
#' @param randRatio2 The randomization ratio between the treatment group and control group for MRCT 2. Defaults to \code{randRatio1}.
#'
#' @details
#' The extended consistency probability via the extended criterion Method I of MHLW (2007) (conditional version),
#' \deqn{\mathrm{Pr}\left(D_{k,\textrm{pool}}\ge \pi D_\textrm{pool}\ |\ T^{(1)}>z_{1-\alpha},T^{(2)}>z_{1-\alpha}\right),}
#' is approximately
#' \deqn{
#'  \frac{1}{(1-\beta_{1})(1-\beta_{2})} \int_{-z_{1-\beta_1}}^{\infty} \int_{-z_{1-\beta_2}}^{\infty} \Phi\left( \frac{\mu_{L|D} (u, v)}{\sigma_{L|D}} \right) \phi(u)\phi(v) \mathrm{d}u\mathrm{d}v,
#' } 
#' where the conditional mean \eqn{\mu_{L|D}} and conditional variance \eqn{\sigma_{L|D}} are given in Qing et al (2025).
#' 
#' Since there is no closed forms of above equations, \code{conProb2} utilizes
#' the \code{\link[stats]{integrate}} function for numerical integration.
#' 
#' Additionally, if we assume equal treatment effects, common regional fractions, common powers and common randomization ratios across 
#' two studies, i.e., \eqn{\sigma^{2(\textrm{t},1)}=\sigma^{2(\textrm{t},2)}, \sigma^{2(\textrm{c},1)}=\sigma^{2(\textrm{c},2)},
#' d^{(1)}=d^{(2)}, f_k^{(1)}=f_k^{(2)}=f_k, \beta_1=\beta_2} and \eqn{r^{(1)}=r^{(2)}}, then the above equation reduces to
#' \deqn{
#'  \frac{1}{(1-\beta)^2} \int_{-z_{1-\beta}}^{\infty} \int_{-z_{1-\beta}}^{\infty} \Phi\left( \frac{(1-\pi) (2z_{1-\alpha} + 2z_{1-\beta} + u + v)}{\sqrt{2(f_k^{-1} - 1)}} \right) \phi(u)\phi(v) \mathrm{d}u \mathrm{d}v.
#' } 
#'  
#' Since now the above equation includes none of \code{d1}=\code{d2}, 
#' \code{sigmaTrt1}=\code{sigmaTrt2}, \code{sigmaCtrl1}=\code{sigmaCtrl2}
#' \code{randRatio1} = \code{randRatio2}, if only the consistency probability
#' is considered, then the values of \code{d1} and \code{sigmaTrt1} could be arbitrary.
#'  
#' The overall sample size is calculated in the same way as \code{\link{conProb}}. 
#'  
#' @returns A list containing the following two components:
#' \describe{
#'   \item{\code{CP}}{The consistency probability, a scalar.}
#'   \item{\code{N1}}{The overall sample size for MRCT 1.}
#'   \item{\code{N2}}{The overall sample size for MRCT 2.}
#' }
#'   
#' @examples 
#' conProb2(alpha = 0.05, power1 = 0.8, power2 = 0.9, rF1 = 0.141, d1 = 1, sigmaTrt1 = 4)
#'
#' 
#' @export
#'
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
  
  # Standard normal quantiles
  za <- qnorm(1 - alpha)
  zb1 <- qnorm(power1)
  zb2 <- qnorm(power2)
  
  # --- Calculate Sample Sizes (Integers) to feed into logic ---
  # Study 1
  num1 <- (sigmaTrt1^2 / randRatio1 + sigmaCtrl1^2) * (za + zb1)^2
  Nc1_int <- ceiling(num1 / d1^2)
  Nt1_int <- ceiling(Nc1_int * randRatio1)
  N1 <- Nt1_int + Nc1_int
  
  # Study 2
  num2 <- (sigmaTrt2^2 / randRatio2 + sigmaCtrl2^2) * (za + zb2)^2
  Nc2_int <- ceiling(num2 / d2^2)
  Nt2_int <- ceiling(Nc2_int * randRatio2)
  N2 <- Nt2_int + Nc2_int
  
  # --- Map Parameters to Logic Variables ---
  beta1 <- 1 - power1
  beta2 <- 1 - power2
  # N1, N2 are set above
  # d1, d2 are set
  # mu_c1, mu_c2
  mu_c1 <- muCtrl1
  mu_c2 <- muCtrl2
  # r1, r2
  r1 <- randRatio1
  r2 <- randRatio2
  # fk1, fk2
  fk1 <- rF1
  fk2 <- rF2
  # sigmas
  sigma_t1 <- sigmaTrt1; sigma_c1 <- sigmaCtrl1
  sigma_t2 <- sigmaTrt2; sigma_c2 <- sigmaCtrl2
  
  # =========================================================
  # LOGIC IMPLEMENTATION (Based on calculate_cp_pool_general)
  # =========================================================
  
  # --- Step 1: Calculate Group Sample Sizes ---
  # Study 1
  N_t1 <- N1 * r1 / (1 + r1)
  N_c1 <- N1 / (1 + r1)
  # Study 2
  N_t2 <- N2 * r2 / (1 + r2)
  N_c2 <- N2 / (1 + r2)
  
  # --- Step 2: Calculate Regional Sub-sample Sizes ---
  Nk_t1 <- fk1 * N_t1; Nk_c1 <- fk1 * N_c1
  Nk_t2 <- fk2 * N_t2; Nk_c2 <- fk2 * N_c2
  
  # --- Step 3: Calculate Means ---
  # mu_t = mu_c + d
  mu_t1 <- mu_c1 + d1
  mu_t2 <- mu_c2 + d2
  
  # --- Step 4: Calculate Weights ---
  # Global pooled weights W^(h,s)
  Denom_T <- N_t1 + N_t2
  Denom_C <- N_c1 + N_c2
  
  W_t1 <- N_t1 / Denom_T; W_t2 <- N_t2 / Denom_T
  W_c1 <- N_c1 / Denom_C; W_c2 <- N_c2 / Denom_C
  
  # Regional pooled weights w_k^(h,s)
  Denom_k_T <- Nk_t1 + Nk_t2
  Denom_k_C <- Nk_c1 + Nk_c2
  
  wk_t1 <- Nk_t1 / Denom_k_T; wk_t2 <- Nk_t2 / Denom_k_T
  wk_c1 <- Nk_c1 / Denom_k_C; wk_c2 <- Nk_c2 / Denom_k_C
  
  # --- Step 5: Calculate Variance Terms ---
  # Var(Y_bar) = sigma^2 / N
  vY_t1 <- sigma_t1^2 / N_t1; vY_c1 <- sigma_c1^2 / N_c1
  vY_t2 <- sigma_t2^2 / N_t2; vY_c2 <- sigma_c2^2 / N_c2
  
  # Variance of effect difference sigma_d^2
  sig2_d1 <- vY_t1 + vY_c1
  sig2_d2 <- vY_t2 + vY_c2
  sig_d1 <- sqrt(sig2_d1)
  sig_d2 <- sqrt(sig2_d2)
  
  # --- Step 6: Calculate Covariance Terms ---
  # lambda_k^(s)
  lam_1 <- wk_t1 * vY_t1 + wk_c1 * vY_c1
  lam_2 <- wk_t2 * vY_t2 + wk_c2 * vY_c2
  
  # Lambda^(s)
  Lam_1 <- W_t1 * vY_t1 + W_c1 * vY_c1
  Lam_2 <- W_t2 * vY_t2 + W_c2 * vY_c2
  
  # C_s = lambda_k - pi * Lambda
  C1 <- lam_1 - pi * Lam_1
  C2 <- lam_2 - pi * Lam_2
  
  # --- Step 7: Calculate E[L] and Var(L) ---
  # E[L] involves means
  Term_Region <- (wk_t1 * mu_t1 - wk_c1 * mu_c1) + (wk_t2 * mu_t2 - wk_c2 * mu_c2)
  Term_Global <- (W_t1 * mu_t1 - W_c1 * mu_c1) + (W_t2 * mu_t2 - W_c2 * mu_c2)
  
  EL <- Term_Region - pi * Term_Global
  
  # Var(L) calculation helper
  calc_var_term <- function(wk, W, vY, fk) {
    (wk^2 * vY / fk) + (pi^2 * W^2 * vY) - (2 * pi * wk * W * vY)
  }
  
  varL <- 0
  varL <- varL + calc_var_term(wk_t1, W_t1, vY_t1, fk1)
  varL <- varL + calc_var_term(wk_c1, W_c1, vY_c1, fk1)
  varL <- varL + calc_var_term(wk_t2, W_t2, vY_t2, fk2)
  varL <- varL + calc_var_term(wk_c2, W_c2, vY_c2, fk2)
  
  # --- Step 8: Calculate Conditional Parameters ---
  # sigma_{L|D}^2 = Var(L) - Sum(C_s^2 / sigma_d^2)
  sig2_LD <- varL - (C1^2 / sig2_d1) - (C2^2 / sig2_d2)
  
  # Check for numerical precision issues
  if (sig2_LD < 0) sig2_LD <- 0
  sig_LD <- sqrt(sig2_LD)
  
  # --- Step 9: Double Numerical Integration ---
  z_b1 <- qnorm(1 - beta1)
  z_b2 <- qnorm(1 - beta2)
  
  # Integrand: Phi( mu_{L|D}(u,v) / sigma_{L|D} ) * phi(u) * phi(v)
  integrand_func <- function(v) {
    sapply(v, function(vi) {
      inner_func <- function(u) {
        # mu_{L|D} = E[L] + (C1/sig_d1)*u + (C2/sig_d2)*v
        mu_val <- EL + (C1 / sig_d1) * u + (C2 / sig_d2) * vi
        pnorm(mu_val / sig_LD) * dnorm(u)
      }
      # Inner integral u: -z_beta1 -> Inf
      res_u <- integrate(inner_func, lower = -z_b1, upper = Inf)$value
      return(res_u * dnorm(vi))
    })
  }
  
  # Outer integral v: -z_beta2 -> Inf
  res_total <- integrate(integrand_func, lower = -z_b2, upper = Inf)$value
  
  # Normalize
  final_cp <- res_total / ((1 - beta1) * (1 - beta2))
  
  return(list(CP = as.numeric(final_cp), N1 = N1, N2 = N2)) 
}