#' Doubly robust estimator for regional average treatment effects
#'
#' @description
#' This function implements a doubly robust estimator for the average treatment 
#' effect (ATE) in a specific region within a MRCT.  
#'
#' @param Y A numeric vector of the observed outcome.
#' @param X A matrix or data frame of covariates.
#' @param D A binary vector (0 or 1) indicating the region (1 for the target region).
#' @param Trt A binary vector (0 or 1) indicating treatment assignment (1 for the treatment group).
#' @param num_boot Integer. The number of bootstrap iterations for inference. Default is 1000.
#' @param conf_level Numeric. Confidence level for the interval. Default is 1000.
#'
#' @return A list containing:
#' \item{estimate}{The point estimate of the regional ATE.}
#' \item{se}{The bootstrap standard error.}
#' \item{conf_int}{A named vector with the lower and upper confidence limits.}
#' \item{p_value_two_sided}{The p-value for a two-sided test (\eqn{H_0: \text{regional ATE} = 0}).}
#' \item{p_value_one_sided}{The p-value for a one-sided test (\eqn{H_0: \text{regional ATE} \leq 0}).}
#'
#' @details
#' The estimator uses the formula:
#' \deqn{\widehat{\tau}_\mathrm{dr} = \frac{1}{n \cdot \mathrm{Pr}(D=1)} \sum_{i=1}^{n} \left[ D_i\{\widehat{m}_1(X_i) - \widehat{m}_0(X_i)\} + \widehat{W}_{1,i}\{Y_i - \widehat{m}_1(X_i)\} - \widehat{W}_{0,i}\{Y_i - \widehat{m}_0(X_i)\} \right]}
#' where \eqn{\widehat{W}_{1,i}} and \eqn{\widehat{W}_{0,i}} are weights derived from the propensity scores of being in the 
#' target region and receiving treatment.
#'
#' @importFrom stats lm predict glm binomial sd qnorm pnorm as.formula
#' @examples  
#' set.seed(42)
#' sampleSize <- 300
#' regionProp <- 0.16
#' sampleSizeRegion <- round(regionProp * sampleSize)
#' sampleSizeNon    <- sampleSize - sampleSizeRegion
#' # define region and treatment (1:1 randomization)
#' region <- c(rep(1, sampleSizeRegion), rep(0, sampleSizeNon))
#' trt <- c(rep(1, sampleSizeRegion/2), rep(0, sampleSizeRegion/2),
#'          rep(1, sampleSizeNon/2),    rep(0, sampleSizeNon/2))
#' # generate covariates with regional heterogeneity
#' x1 <- c(rnorm(sampleSizeRegion, mean = 0.5), rnorm(sampleSizeNon, mean = 1.0))
#' x2 <- rnorm(sampleSize, mean = 0)
#' X  <- data.frame(x1 = x1, x2 = x2)
#' # generate outcome 
#' epsilon <- rnorm(sampleSize, mean = 0)
#' Y <- x1 + x2 + trt * (x1 + x2) + epsilon
#' # estimation regATE_dr
#' results <- regATE_dr(Y = Y, X = X, D = region, Trt = trt, num_boot = 1000)
#' print(results)
#' @export
regATE_dr <- function(Y, X, D, Trt, num_boot = 1000, conf_level = 0.95) {
  
  # --- 1. Data Preparation ---
  # Ensure X is a data frame for formula handling
  df <- data.frame(Y = Y, D = D, Trt = Trt)
  X <- as.data.frame(X)
  cov_names <- colnames(X)
  df <- cbind(df, X)
  
  # Core estimation logic to be called by point estimate and bootstrap
  calc_tau_dr <- function(curr_dat) {
    # Calculate regional proportion (n_region / n_total)
    region_prop <- mean(curr_dat$D)
    
    # A. Outcome Models: E(Y | X, Trt)
    f_y <- stats::as.formula(paste("Y ~", paste(cov_names, collapse = "+")))
    m1 <- stats::lm(f_y, data = curr_dat[curr_dat$Trt == 1, ])
    m0 <- stats::lm(f_y, data = curr_dat[curr_dat$Trt == 0, ])
    
    m1_pred <- stats::predict(m1, newdata = curr_dat)
    m0_pred <- stats::predict(m0, newdata = curr_dat)
    
    # B. Propensity Models: P(D=1 | X) and P(Trt=1 | D, X)
    f_d <- stats::as.formula(paste("D ~", paste(cov_names, collapse = "+")))
    reg_ps_model <- stats::glm(f_d, data = curr_dat, family = stats::binomial)
    reg_ps <- stats::predict(reg_ps_model, newdata = curr_dat, type = "response")
    
    f_t <- stats::as.formula(paste("Trt ~", paste(cov_names, collapse = "+")))
    # Propensity score in target region
    t_ps_reg_mod <- stats::glm(f_t, data = curr_dat[curr_dat$D == 1, ], family = stats::binomial)
    t_ps_reg <- stats::predict(t_ps_reg_mod, newdata = curr_dat, type = "response")
    # Propensity score in other regions
    t_ps_non_mod <- stats::glm(f_t, data = curr_dat[curr_dat$D == 0, ], family = stats::binomial)
    t_ps_non <- stats::predict(t_ps_non_mod, newdata = curr_dat, type = "response")
    
    # C. Weight Construction (r1 = 1, r0 = 1 as per your logic)
    resid1 <- curr_dat$Y - m1_pred
    resid0 <- curr_dat$Y - m0_pred
    
    w1_num <- curr_dat$D * curr_dat$Trt * reg_ps + (1 - curr_dat$D) * curr_dat$Trt * reg_ps
    w1_den <- reg_ps * t_ps_reg + (1 - reg_ps) * t_ps_non
    W1 <- w1_num / w1_den
    
    w0_num <- curr_dat$D * (1 - curr_dat$Trt) * reg_ps + (1 - curr_dat$D) * (1 - curr_dat$Trt) * reg_ps
    w0_den <- reg_ps * (1 - t_ps_reg) + (1 - reg_ps) * (1 - t_ps_non)
    W0 <- w0_num / w0_den
    
    # D. DR Estimator Calculation
    dr_val <- (1 / region_prop) * mean(curr_dat$D * (m1_pred - m0_pred) + W1 * resid1 - W0 * resid0)
    return(dr_val)
  }
  
  # --- 2. Point Estimate ---
  point_est <- calc_tau_dr(df)
  
  # --- 3. Stratified Bootstrap for Inference ---
  boot_dist <- numeric(num_boot)
  idx_1 <- which(df$D == 1)
  idx_0 <- which(df$D == 0)
  
  for (b in 1:num_boot) {
    # Stratify by D to ensure region proportion remains constant
    b_idx <- c(
      sample(idx_1, length(idx_1), replace = TRUE),
      sample(idx_0, length(idx_0), replace = TRUE)
    )
    boot_dist[b] <- suppressWarnings(calc_tau_dr(df[b_idx, ]))
  }
  
  # --- 4. Statistical Inference ---
  se <- stats::sd(boot_dist, na.rm = TRUE)
  z_stat <- point_est / se
  
  alpha <- 1 - conf_level
  z_crit <- stats::qnorm(1 - alpha / 2)
  ci <- c(lower = point_est - z_crit * se, upper = point_est + z_crit * se)
  
  p_two <- 2 * (1 - stats::pnorm(abs(z_stat)))
  p_one <- 1 - stats::pnorm(z_stat) # H0: tau <= 0
  
  return(list(
    estimate = point_est,
    se = se,
    conf_int = ci,
    p_value_two_sided = p_two,
    p_value_one_sided = p_one
  ))
}
 