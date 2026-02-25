#' @name conProbII
#'
#' @title Consistency probability for MRCT via criterion Method II of MHLW (2007) (conditional version) under fixed effects model
#'
#' @description Calculate the consistency probability for MRCT via criterion Method II of MHLW (2007) (conditional version) under fixed effects model.
#'
#' @param alpha The Type I error.
#' @param power Power. 
#' @param rF The regional fractions, a \code{K x 1} vector with each component 
#' representing the regional fraction in region \eqn{k} and the sum must equal 1.
#' @param d The true mean of difference of response.
#' @param sigmaTrt The standard deviation of response in the treatment group.
#' @param sigmaCtrl The standard deviation of response in the control group. Defaults to \code{sigmaTrt}.
#' @param pTrt The mean of the response in the treatment group.
#' @param pCtrl The mean of the response in the control group.
#' @param randRatio The randomization ratio between the treatment group and control group. Defaults to 1.
#' @param responseType The type of response. One of "continuous" and "binary".
#' @param B The number of simulation by Monto Carlo for \code{responseType} = "binary". Defaults to 100,000.
#'
#' @details
#' The consistency probability via criterion Method II of MHLW (2007) (conditional version),
#' \deqn{\mathrm{Pr}\left(D_{k}\geq 0, k=1,\ldots,K\ |\ T>z_{1-\alpha}\right),}
#' is approximately
#' \deqn{
#' \frac{1}{1-\beta}\int_{-z_{1-\beta}}^{\infty} \Phi_K\left\{ \mathbf{0};-(u + z_{1-\alpha} + z_{1-\beta})\mathbf{1}_K, \mathbf{\Lambda} - \mathbf{1}_K \mathbf{1}_K^\top \right\} \phi(u) \mathrm{d} u,
#' }
#' where \eqn{\Phi_K(\mathbf{x}; \vec{\mu}, \Sigma)} denotes the cumulative distribution function of a \eqn{K}-dimensional normal distribution with mean \eqn{\vec{\mu}}
#' and covariance matrix \eqn{\Sigma}, \eqn{\mathbf{1}_K} is the vector of \eqn{K} ones, \eqn{\mathbf{\Lambda} = \text{diag}(f_{1}^{-1}, \dots, f_{K}^{-1})}. 
#' 
#' Since there is no closed form of above equation, \code{conProbII} utilizes the \code{\link[mvtnorm]{pmvnorm}} function for numerical integration.
#'   
#' As the first equation includes none of \code{d}, \code{sigmaTrt}, 
#' \code{sigmaCtrl} and \code{randRatio}, if only the consistency probability
#' is considered, then the values of \code{d} and \code{sigmaTrt} could be arbitrary.
#' 
#' For binary response, the above approximation loses precision under moderate sample size.
#' Hence \code{conProbII} applies Monto Carlo to calculate the correct consistency probability.
#'
#' The overall sample size is calculated in the same way as \code{\link{conProb}}. 
#' 
#' @returns A list containing the following two components:
#' \describe{
#'   \item{\code{CP}}{The consistency probability, a scalar.}
#'   \item{\code{N}}{The overall sample size.}
#' } 
#' 
#' @export 
#' @importFrom mvtnorm pmvnorm GenzBretz
#' @examples
#' # Example 1 in Qing et al (2025)
#' alpha <- 0.05
#' power <- 0.8
#' conProbII(alpha, power, rF = rep(1/2,2), d = 1, sigmaTrt = 4, responseType = "continuous")
#' conProbII(alpha, power, rF = rep(1/3,3), d = 1, sigmaTrt = 4, responseType = "continuous")
#' conProbII(alpha, power, rF = rep(1/4,4), d = 1, sigmaTrt = 4, responseType = "continuous")
#' rFk1 <- 0.110
#' rFk2 <- rFk3 <- (1-rFk1)/2
#' rF <- c(rFk1, rFk2, rFk3)
#' conProbII(alpha, power, rF = rF, d = 1, sigmaTrt = 4, responseType = "continuous") 
#' # Example 2 in Qing et al (2025)
#' rFk1 <- 0.149
#' rFk2 <- rFk3 <- (1-rFk1)/2
#' rF <- c(rFk1, rFk2, rFk3)
#' set.seed(123)
#' conProbII(alpha, power, rF = rF, pTrt = 0.8, pCtrl = 0.7, responseType = "binary")
#' rFk1 <- 0.14
#' rFk2 <- rFk3 <- (1-rFk1)/2
#' rF <- c(rFk1, rFk2, rFk3)
#' set.seed(123)
#' conProbII(alpha, power, rF = rF, pTrt = 0.7, pCtrl = 0.6, responseType = "binary")
#' rFk1 <- 0.110
#' rFk2 <- rFk3 <- (1-rFk1)/2
#' rF <- c(rFk1, rFk2, rFk3)
#' set.seed(123)
#' conProbII(alpha, power, rF = rF, pTrt = 0.8, pCtrl = 0.7, responseType = "binary")
#' conProbII(alpha, power, rF = rF, pTrt = 0.7, pCtrl = 0.6, responseType = "binary")
#'
#' @export
#'
conProbII <- function(alpha, 
                    power,  
                    rF, 
                    d,  
                    sigmaTrt, 
                    sigmaCtrl = sigmaTrt,  
                    pTrt,
                    pCtrl,
                    randRatio = 1,
                    responseType = c("continuous","binary"),
                    B = 100000){
  if ((sum(rF) != 1) || any(rF <= 0) ){
    stop("the sum of 'rF' must equal 1 and each component must be greater than 0")
    } 
  ### consistency probability
  if (responseType == "continuous"){
    z_alpha = qnorm(1 - alpha)
    z_beta = qnorm(power)
    K <- length(rF)
    Lambda = diag(1 / rF)
    Sigma = Lambda - matrix(1, K, K)
    
    integrand <- function(z) {
      sapply(z, function(zi) {
        mu_val = -(zi + z_alpha + z_beta)
        mu_vec = rep(mu_val, K) 
        mvtnorm::pmvnorm(upper = -mu_vec, mean = rep(0, K), sigma = Sigma, algorithm = GenzBretz(abseps = 1e-5))[1] * dnorm(zi)
      })
    }
    CP <- integrate(integrand, lower = -z_beta, upper = Inf, rel.tol = 1e-4)$value/(power)
  ### sample size
  NTrt <-  ceiling((1/randRatio*sigmaTrt^2 + sigmaCtrl^2)*(qnorm(1-alpha) + qnorm(power))^2/d^2)
  NCtrl <- ceiling(randRatio*NTrt)
  N <- NTrt + NCtrl
  } else if (responseType == "binary"){
    K <- length(rF)
    d <- pTrt - pCtrl
    sigmaTrt <- sqrt(pTrt*(1-pTrt))
    sigmaCtrl <- sqrt(pCtrl*(1-pCtrl))
    NTrt <-  ceiling((1/randRatio*sigmaTrt^2 + sigmaCtrl^2)*(qnorm(1-alpha) + qnorm(power))^2/d^2)
    NCtrl <- ceiling(randRatio*NTrt)
    NTrtReg <- ceiling(rF*NTrt)
    NCtrlReg <- ceiling(rF*NCtrl)
    NTrt <- sum(NTrtReg)
    NCtrl <- sum(NCtrlReg)
    N <- NTrt + NCtrl
    consistencyCount <- 0
    for (b in 1:B){
      yTrt <- rbinom(NTrt,size = 1, prob = pTrt) 
      yCtrl <- rbinom(NCtrl,size = 1, prob = pCtrl)  
      pTrtHat <- sum(yTrt)/NTrt
      pCtrlHat <- sum(yCtrl)/NCtrl
      sigmaDHat <- sqrt((randRatio + 1)*(pTrtHat*(1-pTrtHat)+randRatio*pCtrlHat*(1-pCtrlHat))/randRatio/N)
      dHat <- pTrtHat - pCtrlHat 
      pTrtRegHat <- rep(0,K)
      pCtrlRegHat <- rep(0,K)
      startRowTrt <- startRowCtrl <- 1
      for (k in 1:K){
        endRowTrt <- NTrtReg[k] + startRowTrt - 1
        endRowCtrl <- NCtrlReg[k] + startRowCtrl - 1 
        pTrtRegHat[k] <- sum(yTrt[startRowTrt:endRowTrt])/NTrtReg[k]
        pCtrlRegHat[k] <- sum(yCtrl[startRowCtrl:endRowCtrl])/NCtrlReg[k]
        startRowTrt <- endRowTrt + 1
        startRowCtrl <- endRowCtrl + 1
      }
      dRegHat <- pTrtRegHat - pCtrlRegHat
      if ( all(dRegHat > 0) && (dHat/sigmaDHat > qnorm(1-alpha)) ){
        consistencyCount <- consistencyCount + 1
      }
    }
    CP <- as.numeric(consistencyCount/B)/power 
  } else {stop("'responseType' must be one of 'continuous' and 'binary'")} 
  return(list(CP = CP, N = N))
}
