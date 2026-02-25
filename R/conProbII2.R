#' @name conProbII2
#'
#' @title Consistency probability for two MRCTs via criterion Method II (conditional version)
#'
#' @description Calculate the consistency probability for two MRCTs via criterion Method II of MHLW (2007) (conditional version) under fixed effects model.
#'
#' @param alpha The Type I error.
#' @param power1 Power for MRCT 1.
#' @param power2 Power for MRCT 2. Defaults to \code{power1}.
#' @param rF1 The regional fractions for MRCT 1, a \code{K x 1} vector with each component 
#' representing the regional fraction in region \eqn{k} and the sum must equal 1.
#' @param rF2 The regional fractions for MRCT 2, a \code{K x 1} vector with each component 
#' representing the regional fraction in region \eqn{k} and the sum must equal 1. Defaults to \code{rF1}.
#' @param d1 The true mean of difference of response for MRCT 1.
#' @param d2 The true mean of difference of response for MRCT 2. Defaults to \code{d1}.
#' @param sigmaTrt1 The standard deviation of response in the treatment group for MRCT 1.
#' @param sigmaCtrl1 The standard deviation of response in the control group for MRCT 1. Defaults to \code{sigmaTrt1}.
#' @param sigmaTrt2 The standard deviation of response in the treatment group for MRCT 2. Defaults to \code{sigmaTrt1}.
#' @param sigmaCtrl2 The standard deviation of response in the control group for MRCT 2. Defaults to \code{sigmaTrt2}.
#' @param pTrt1 The mean of the response in the treatment group for MRCT 1.
#' @param pCtrl1 The mean of the response in the control group for MRCT 1.
#' @param pTrt2 The mean of the response in the treatment group for MRCT 2. Defaults to \code{pTrt1}.
#' @param pCtrl2 The mean of the response in the control group for MRCT 2. Defaults to \code{pTrt2}.
#' @param randRatio1 The randomization ratio between the treatment group and control group for MRCT 1.
#' @param randRatio2 The randomization ratio between the treatment group and control group for MRCT 2. Defaults to \code{randRatio1}.
#' @param responseType The type of response. One of "continuous" and "binary".
#' @param B The number of simulation by Monto Carlo for \code{responseType} = "binary". Defaults to 100,000.
#' 
#' @details
#' The extended consistency probability via the extended riterion Method II of MHLW (2007) (conditional version),
#' \deqn{\mathrm{Pr}\left(D_{k,\textrm{pool}}\geq 0, k=1,\cdots,K\ |\ T^{(1)}>z_{1-\alpha},T^{(2)}>z_{1-\alpha}\right)>1-\gamma,}
#' is approximately
#' \deqn{
#' \frac{1}{(1-\beta_{1})(1-\beta_{2})}\int_{-z_{1-\beta_{1}}}^{\infty}\int_{-z_{1-\beta_{2}}}^{\infty} \Phi_K\left( \mathbf{0}; -\vec{\xi}(u,v), \Omega \right) \phi(u)\phi(v) \mathrm{d}u \mathrm{d}v,
#' }
#' where \eqn{\vec{\xi}} and \eqn{\Omega} are given in Qing et al (2025).
#' 
#' Since there is no closed form of above equation, \code{conProbII2} utilizes the \code{\link[mvtnorm]{pmvnorm}} function for numerical integration.
#' 
#' For binary response, the above approximation loses precision under moderate sample size.
#' Hence \code{conProbII2} applies Monto Carlo to calculate the correct consistency probability.
#'
#' The overall sample size is calculated in the same way as \code{\link{conProb}}.
#' 
#' @returns A list containing the following two components:
#' \describe{
#'   \item{\code{CP}}{The consistency probability, a scalar.}
#'   \item{\code{N1}}{The overall sample size for MRCT 1.}
#'   \item{\code{N2}}{The overall sample size for MRCT 2.}
#' }
#' @examples
#' # Example 3 in Qing et al (2025)
#' alpha <- 0.05
#' power1 <- 0.8
#' conProbII2(alpha = alpha, power1 = power1, rF1 = rep(1/2,2), d1 = 1, 
#' sigmaTrt1 = 4, responseType = "continuous")
#' conProbII2(alpha = alpha, power1 = power1, rF1 = rep(1/3,3), d1 = 1, 
#' sigmaTrt1 = 4, responseType = "continuous")
#' conProbII2(alpha = alpha, power1 = power1, rF1 = rep(1/4,4), d1 = 1, 
#' sigmaTrt1 = 4, responseType = "continuous")
#' rFk1 <- 0.044
#' rFk2 <- rFk3 <- (1-rFk1)/2
#' rF <- c(rFk1, rFk2, rFk3)
#' conProbII2(alpha = alpha, power1 = power1, rF1 = rF, d1 = 1, 
#' sigmaTrt1 = 4, responseType = "continuous")
#' 
#' # Example 4 in Qing et al (2025)
#' alpha <- 0.05
#' power1 <- 0.8
#' rFk1 <- 0.06
#' rFk2 <- rFk3 <- (1-rFk1)/2
#' rF <- c(rFk1, rFk2, rFk3)
#' conProbII2(alpha = alpha, power1 = power1, rF1 = rF, pTrt1 = 0.9, 
#' pCtrl1 = 0.8, responseType = "binary")
#' rFk1 <- 0.044
#' rFk2 <- rFk3 <- (1-rFk1)/2
#' rF <- c(rFk1, rFk2, rFk3)
#' conProbII2(alpha = alpha, power1 = power1, rF1 = rF, pTrt1 = 0.9, 
#' pCtrl1 = 0.8, responseType = "binary")
#'
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
                       B = 100000){
  
  if ((sum(rF1) != 1) || any(rF1 <= 0) ){
    stop("the sum of 'rF1' must equal 1 and each component must be greater than 0")
  } 
  if ((sum(rF2) != 1) || any(rF2 <= 0) ){
    stop("the sum of 'rF2' must equal 1 and each component must be greater than 0")
  } 
  if ( length(rF1) != length(rF2)){
    stop("the length of 'rF1' must equal the length of 'rF2'")
  }
  
  # Check for required package
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Package 'mvtnorm' is required. Please install it.")
  }
  
  ### consistency probability
  if (responseType == "continuous"){
    # 1. Total Sample Sizes (N1, N2)
    # N = (r+1) * (sigma_t^2/r + sigma_c^2) * (z_alpha + z_beta)^2 / d^2
    N1 <- (randRatio1+1)*(sigmaTrt1^2/randRatio1+sigmaCtrl1^2)*(qnorm(1-alpha)+qnorm(power1))^2/(d1^2)
    N2 <- (randRatio2+1)*(sigmaTrt2^2/randRatio2+sigmaCtrl2^2)*(qnorm(1-alpha)+qnorm(power2))^2/(d2^2)
    
    # 2. Derived Sample Size Parameters
    Nt1 <- N1 * randRatio1 / (randRatio1 + 1)
    Nc1 <- N1 * 1 / (randRatio1 + 1)
    Nt2 <- N2 * randRatio2 / (randRatio2 + 1)
    Nc2 <- N2 * 1 / (randRatio2 + 1)
    
    # Variances of Study Means (Var(Y_bar))
    var_Ybar_t1 <- sigmaTrt1^2 / Nt1
    var_Ybar_c1 <- sigmaCtrl1^2 / Nc1
    var_Ybar_t2 <- sigmaTrt2^2 / Nt2
    var_Ybar_c2 <- sigmaCtrl2^2 / Nc2
    
    # Variance of Treatment Effect Estimate D^(s) (sigma_d^2)
    sigmad1_sq <- var_Ybar_t1 + var_Ybar_c1
    sigmad2_sq <- var_Ybar_t2 + var_Ybar_c2
    
    # 3. Regional Weights and Lambda Calculation
    K <- length(rF1)
    
    # Regional sample sizes (Treatment and Control)
    Nt1_k <- rF1 * Nt1
    Nc1_k <- rF1 * Nc1
    Nt2_k <- rF2 * Nt2
    Nc2_k <- rF2 * Nc2
    
    # Pooling weights for each region k: w_k^(h,s)
    # w_k^(t,1) = N_k^(t,1) / (N_k^(t,1) + N_k^(t,2))
    denom_t <- Nt1_k + Nt2_k
    denom_c <- Nc1_k + Nc2_k
    
    w_t1 <- Nt1_k / denom_t
    w_t2 <- Nt2_k / denom_t
    w_c1 <- Nc1_k / denom_c
    w_c2 <- Nc2_k / denom_c
    
    # Lambda_k^(s) = Cov(D_{k,pool}, D^(s))
    lambda1 <- w_t1 * var_Ybar_t1 + w_c1 * var_Ybar_c1
    lambda2 <- w_t2 * var_Ybar_t2 + w_c2 * var_Ybar_c2
    
    # 4. Joint Covariance Matrix construction
    # Vector X = (D_{1,pool}, ..., D_{K,pool}, Z^(1), Z^(2))
    # Dimensions: K + 2
    
    # Block 1: Cov(D_pool, D_pool)
    # Independent regions -> Diagonal matrix
    var_Dk_pool <- (w_t1^2 * sigmaTrt1^2 / Nt1_k) + (w_c1^2 * sigmaCtrl1^2 / Nc1_k) +
      (w_t2^2 * sigmaTrt2^2 / Nt2_k) + (w_c2^2 * sigmaCtrl2^2 / Nc2_k)
    Sigma_DD <- diag(var_Dk_pool, nrow=K, ncol=K)
    
    # Block 2: Cov(Z, Z) -> Identity 2x2
    Sigma_ZZ <- diag(1, 2, 2)
    
    # Block 3: Cov(D_pool, Z)
    # Cov(D_{k,pool}, Z^{(s)}) = Cov(D_{k,pool}, D^{(s)}) / sigma_d^{(s)} = lambda_k^{(s)} / sigma_d^{(s)}
    Cov_DZ1 <- lambda1 / sqrt(sigmad1_sq) # Vector length K
    Cov_DZ2 <- lambda2 / sqrt(sigmad2_sq) # Vector length K
    
    # Assemble Sigma Total (K+2 x K+2)
    # [ Sigma_DD   Cov_DZ1  Cov_DZ2 ]
    # [ Cov_DZ1^T     1        0    ]
    # [ Cov_DZ2^T     0        1    ]
    
    Sigma_total <- matrix(0, nrow = K+2, ncol = K+2)
    Sigma_total[1:K, 1:K] <- Sigma_DD
    Sigma_total[(K+1), (K+1)] <- 1
    Sigma_total[(K+2), (K+2)] <- 1
    
    # Fill off-diagonals
    Sigma_total[1:K, K+1] <- Cov_DZ1
    Sigma_total[K+1, 1:K] <- Cov_DZ1
    Sigma_total[1:K, K+2] <- Cov_DZ2
    Sigma_total[K+2, 1:K] <- Cov_DZ2
    
    # 5. Mean Vector
    # Mean(D_{k,pool}) = w_t1*d1 + w_t2*d2 (Assuming baselines cancel or w_t=w_c)
    mean_D <- w_t1 * d1 + w_t2 * d2
    # Mean(Z) = 0 (Standardized)
    Mean_total <- c(mean_D, 0, 0)
    
    # 6. Integration Limits
    # We want P(D_k >= 0, Z^(1) > -z_(1-beta1), Z^(2) > -z_(1-beta2))
    lower_limits <- c(rep(0, K), -qnorm(power1), -qnorm(power2))
    upper_limits <- rep(Inf, K+2)
    
    # 7. Calculate Probability
    # Use pmvnorm directly on the joint distribution
    prob_joint <- mvtnorm::pmvnorm(lower = lower_limits, 
                                   upper = upper_limits, 
                                   mean = Mean_total, 
                                   sigma = Sigma_total)
    
    # Conditional Probability: Joint / Denominator
    # Denominator approx power1 * power2
    CP <- as.numeric(prob_joint) / (power1 * power2)
    
    ### Sample size return variables (integers)
    NTrt1 <-  ceiling((1/randRatio1*sigmaTrt1^2 + sigmaCtrl1^2)*(qnorm(1-alpha) + qnorm(power1))^2/d1^2)
    NCtrl1 <- ceiling(randRatio1*NTrt1)
    N1_out <- NTrt1 + NCtrl1
    NTrt2 <-  ceiling((1/randRatio2*sigmaTrt2^2 + sigmaCtrl2^2)*(qnorm(1-alpha) + qnorm(power2))^2/d2^2)
    NCtrl2 <- ceiling(randRatio2*NTrt2)
    N2_out <- NTrt2 + NCtrl2
    
    return(list(CP = CP, N1 = N1_out, N2 = N2_out))
    
  } else if (responseType == "binary"){
    # Binary case remains simulation-based
    K <- length(rF1)
    d1 <- pTrt1 - pCtrl1
    d2 <- pTrt2 - pCtrl2
    sigmaTrt1 <- sqrt(pTrt1*(1-pTrt1))
    sigmaCtrl1 <- sqrt(pCtrl1*(1-pCtrl1))
    sigmaTrt2 <- sqrt(pTrt2*(1-pTrt2))
    sigmaCtrl2 <- sqrt(pCtrl2*(1-pCtrl2))
    #
    NTrt1 <-  ceiling((1/randRatio1*sigmaTrt1^2 + sigmaCtrl1^2)*(qnorm(1-alpha) + qnorm(power1))^2/d1^2)
    NCtrl1 <- ceiling(randRatio1*NTrt1)
    NTrtReg1 <- ceiling(rF1*NTrt1)
    NCtrlReg1 <- ceiling(rF1*NCtrl1)
    NTrt1 <- sum(NTrtReg1)
    NCtrl1 <- sum(NCtrlReg1)
    N1 <- NTrt1 + NCtrl1
    #
    NTrt2 <-  ceiling((1/randRatio2*sigmaTrt2^2 + sigmaCtrl2^2)*(qnorm(1-alpha) + qnorm(power2))^2/d2^2)
    NCtrl2 <- ceiling(randRatio2*NTrt2)
    NTrtReg2 <- ceiling(rF2*NTrt2)
    NCtrlReg2 <- ceiling(rF2*NCtrl2)
    NTrt2 <- sum(NTrtReg2)
    NCtrl2 <- sum(NCtrlReg2)
    N2 <- NTrt2 + NCtrl2
    consistencyCount <- 0
    
    # Binary simulation loop
    for (b in 1:B){
      yTrt1 <- rbinom(NTrt1,size = 1, prob = pTrt1) 
      yCtrl1 <- rbinom(NCtrl1,size = 1, prob = pCtrl1)  
      pTrtHat1 <- sum(yTrt1)/NTrt1
      pCtrlHat1 <- sum(yCtrl1)/NCtrl1
      sigmaDHat1 <- sqrt((randRatio1 + 1)*(pTrtHat1*(1-pTrtHat1)+randRatio1*pCtrlHat1*(1-pCtrlHat1))/randRatio1/N1)
      dHat1 <- pTrtHat1 - pCtrlHat1 
      numTrtRegHat1 <- rep(0,K)
      numCtrlRegHat1 <- rep(0,K)
      startRowTrt1 <- startRowCtrl1 <- 1
      for (k in 1:K){
        endRowTrt1 <- NTrtReg1[k] + startRowTrt1 - 1
        endRowCtrl1 <- NCtrlReg1[k] + startRowCtrl1 - 1 
        numTrtRegHat1[k] <- sum(yTrt1[startRowTrt1:endRowTrt1]) 
        numCtrlRegHat1[k] <- sum(yCtrl1[startRowCtrl1:endRowCtrl1]) 
        startRowTrt1 <- endRowTrt1 + 1
        startRowCtrl1 <- endRowCtrl1 + 1
      }
      #
      yTrt2 <- rbinom(NTrt2,size = 1, prob = pTrt2) 
      yCtrl2 <- rbinom(NCtrl2,size = 1, prob = pCtrl2)  
      pTrtHat2 <- sum(yTrt2)/NTrt2
      pCtrlHat2 <- sum(yCtrl2)/NCtrl2
      sigmaDHat2 <- sqrt((randRatio2 + 1)*(pTrtHat2*(1-pTrtHat2)+randRatio2*pCtrlHat2*(1-pCtrlHat2))/randRatio2/N2)
      dHat2 <- pTrtHat2 - pCtrlHat2 
      numTrtRegHat2 <- rep(0,K)
      numCtrlRegHat2 <- rep(0,K)
      startRowTrt2 <- startRowCtrl2 <- 1
      for (k in 1:K){
        endRowTrt2 <- NTrtReg2[k] + startRowTrt2 - 1
        endRowCtrl2 <- NCtrlReg2[k] + startRowCtrl2 - 1 
        numTrtRegHat2[k] <- sum(yTrt2[startRowTrt2:endRowTrt2]) 
        numCtrlRegHat2[k] <- sum(yCtrl2[startRowCtrl2:endRowCtrl2]) 
        startRowTrt2 <- endRowTrt2 + 1
        startRowCtrl2 <- endRowCtrl2 + 1
      }
      numTrtRegHat <- numTrtRegHat1 + numTrtRegHat2
      numCtrlRegHat <- numCtrlRegHat1 + numCtrlRegHat2
      dRegHat <- numTrtRegHat/(NTrtReg1 + NTrtReg2) - numCtrlRegHat/(NCtrlReg1 + NCtrlReg2)
      if ( all(dRegHat > 0) && (dHat1/sigmaDHat1 > qnorm(1-alpha)) && (dHat2/sigmaDHat2 > qnorm(1-alpha))){
        consistencyCount <- consistencyCount + 1
      }
    }
    CP <- as.numeric(consistencyCount/B)/power1/power2 
    return(list(CP = CP, N1 = N1, N2 = N2))
  } else {
    stop("'responseType' must be one of 'continuous' and 'binary'")
  } 
}
