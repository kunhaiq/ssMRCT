#' @name conProbII2
#'
#' @title Consistency probability for two MRCTs via Japan's criterion II (conditional version)
#'
#' @description Calculate the consistency probability for two MRCTs via Japan's criterion II (conditional version).
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
#' The extended consistency probability via the extended Japan's criterion II (conditional version),
#' \deqn{\mathrm{Pr}\left(D_{k,\textrm{pool}}\geq 0, k=1,\cdots,K\ |\ T^{(1)}>z_{1-\alpha},T^{(2)}>z_{1-\alpha}\right),}
#' is approximately
#' \deqn{
#'  \begin{aligned} 
#'    &\frac{1}{(1-\beta_{1})(1-\beta_{2})}\int_{-z_{1-\beta_{1}}}^{\infty}\int_{-z_{1-\beta_{2}}}^{\infty}  \notag\\ 
#'    &\quad\prod_{k=1}^{K}\Phi\left(\frac{ w^{(1)}\sigma^{(1)}_{d}u+w^{(2)}\sigma^{(2)}_{d}v + w^{(1)}d^{(1)} + w^{(2)}d^{(2)} }{\sqrt{\left\{\left(f_{k}^{(1)}\right)^{-1}-1\right\} \left(w^{(1)}\sigma^{(1)}_{d}\right)^{2} + \left\{\left(f_{k}^{(2)}\right)^{-1}-1\right\} \left(w^{(2)}\sigma^{(2)}_{d}\right)^{2}}} \right)\phi(u)\phi(v)dudv,
#' \end{aligned} 
#' }
#' where \eqn{w^{(s)}} and \eqn{\sigma_d^{2}} are defined in \code{\link{conProb2}}.
#'  
#' Since there is no closed forms of above equations, \code{conProbII2} utilizes
#' the \code{\link[cubature]{adaptIntegrate}} function for numerical integration.
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
#'   
#' @examples
#'
#' ### Example 3
#' alpha <- 0.05
#' power1 <- 0.8 
#' conProbII2(alpha, power1, rF1 = rep(1/2,2), d1 = 1, sigmaTrt1 = 4, responseType = "continuous")
#' conProbII2(alpha, power1, rF1 = rep(1/3,3), d1 = 1, sigmaTrt1 = 4, responseType = "continuous")
#' conProbII2(alpha, power1, rF1 = rep(1/4,4), d1 = 1, sigmaTrt1 = 4, responseType = "continuous") 
#' rF11 <- 0.044
#' rF12 <- rF13 <- (1-rF11)/2
#' rF1 <- c(rF11, rF12, rF13)
#' conProbII2(alpha, power1, rF1 = rF1, d1 = 1, sigmaTrt1 = 4, responseType = "continuous") 
#' 
#' ### Example 4
#' rF11 <- 0.06
#' rF12 <- rF13 <- (1-rF11)/2
#' rF1 <- c(rF11, rF12, rF13)
#' set.seed(123)
#' conProbII2(alpha, power1, rF1 = rF1, pTrt1 = 0.8, pCtrl1 = 0.7, responseType = "binary") 
#' rF11 <- 0.044
#' rF12 <- rF13 <- (1-rF11)/2
#' rF1 <- c(rF11, rF12, rF13)
#' set.seed(123)
#' conProbII2(alpha, power1, rF1 = rF1, pTrt1 = 0.8, pCtrl1 = 0.7, responseType = "binary") 
#'
#' @export
#'
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
  ### consistency probability
  if (responseType == "continuous"){
    N1 <- (randRatio1+1)*(sigmaTrt1^2/randRatio1+sigmaCtrl1^2)*(qnorm(1-alpha)+qnorm(power1))^2/(d1^2)
    N2 <- (randRatio2+1)*(sigmaTrt2^2/randRatio2+sigmaCtrl2^2)*(qnorm(1-alpha)+qnorm(power2))^2/(d2^2)
    sigmad1 <- sqrt((randRatio1+1)*(sigmaTrt1^2/randRatio1+sigmaCtrl1^2)/N1)
    sigmad2 <- sqrt((randRatio2+1)*(sigmaTrt2^2/randRatio2+sigmaCtrl2^2)/N2)
    w1 <- N1/(N1+N2)
    w2 <- N2/(N1+N2)
    rF <- as.matrix(data.frame(rF1 = rF1, rF2 = rF2))
    integerandElement <- function(x,rF1k,rF2k){
      return(pnorm((w1*sigmad1*x[1]+w2*sigmad2*x[2]+w1*d1+w2*d2)/sqrt((1/rF1k-1)*w1^2*sigmad1^2+(1/rF2k-1)*w2^2*sigmad2^2)))
    } 
    integerand <- function(x){
      resultsVector <- apply(rF, MARGIN = 1, FUN = function(rF_row) {
        integerandElement(x, rF1k = rF_row[1], rF2k = rF_row[2])
      })
      value <- as.numeric(prod(resultsVector)*dnorm(x[1])*dnorm(x[2]))
      return(value)
    } 
    integerandVec <- Vectorize(integerand)
    CP <- cubature::adaptIntegrate(integerand,lower=c(-qnorm(power1),-qnorm(power2)),upper=c(Inf,Inf))
    CP <- CP$integral/(power1*power2)
    CP <- as.numeric(CP)
    ### sample size 
    NTrt1 <-  ceiling((1/randRatio1*sigmaTrt1^2 + sigmaCtrl1^2)*(qnorm(1-alpha) + qnorm(power1))^2/d1^2)
    NCtrl1 <- ceiling(randRatio1*NTrt1)
    N1 <- NTrt1 + NCtrl1
    NTrt2 <-  ceiling((1/randRatio2*sigmaTrt2^2 + sigmaCtrl2^2)*(qnorm(1-alpha) + qnorm(power2))^2/d2^2)
    NCtrl2 <- ceiling(randRatio2*NTrt2)
    N2 <- NTrt2 + NCtrl2
  } else if (responseType == "binary"){
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
  } else {stop("'responseType' must be one of 'continuous' and 'binary'")} 
  return(list(CP = CP, N1 = N1, N2 = N2))
}
