#' @name conProb2MRCTs
#'
#' @title Consistency probability for two MRCTs in general
#'
#' @description Calculate the (conditional) consistency probability for two MRCTs in general
#'
#' @param alpha The type I error.
#' @param power1 Power for MRCT 1.
#' @param power2 Power for MRCT 2.
#' @param pi The threshold ratio in Japan's criterion I (conditional version).
#' @param rF1 The regional fraction for MRCT 1.
#' @param rF2 The regional fraction for MRCT 2.
#' @param d1 The true mean of difference of response for MRCT 1.
#' @param d2 The true mean of difference of response for MRCT 2.
#' @param sigmaTrt1 The standard deviation of response in the treatment group for MRCT 1.
#' @param sigmaCtrl1 The standard deviation of response in the control group for MRCT 1.
#' @param sigmaTrt2 The standard deviation of response in the treatment group for MRCT 2.
#' @param sigmaCtrl2 The standard deviation of response in the control group for MRCT 2.
#' @param randRatio1 The randomization ratio between the treatment group and control group for MRCT 1.
#' @param randRatio2 The randomization ratio between the treatment group and control group for MRCT 2.
#'
#' @details
#' The extended (conditional) consistency probability, \eqn{\mathrm{Pr}\left(D_{k,\textrm{pool}}\ge \pi D_\textrm{pool}\ |\ T^{(1)}>z_{1-\alpha},T^{(2)}>z_{1-\alpha}\right)},
#' equals
#' \deqn{
#'  \begin{aligned}
#'  &\frac{1}{(1-\beta_{1})(1-\beta_{2})}\int_{-z_{1-\beta_{1}}}^{\infty}\int_{-z_{1-\beta_{2}}}^{\infty} \\
#'  &\quad
#'  \Phi\left(\frac{(1-\pi)\left\{w^{(1)}\sigma^{(1)}_{d}u+w^{(2)}\sigma^{(2)}_{d}v + w^{(1)}d^{(1)} + w^{(2)}d^{(2)}\right\}}{\sqrt{\left\{\left(f_{k}^{(1)}\right)^{-1}-1\right\} \left(w^{(1)}\sigma^{(1)}_{d}\right)^{2} + \left\{\left(f_{k}^{(2)}\right)^{-1}-1\right\} \left(w^{(2)}\sigma^{(2)}_{d}\right)^{2}}} \right)\phi(u)\phi(v)dudv.
#'  \end{aligned}
#' }
#' Since there is no closed form of above equation, \code{conProb2MRCTs} utilizes
#' the \code{\link[cubature]{adaptIntegrate}} function for numerical integration.
#'
#' @examples
#' \dontrun{
#' conProb2MRCTs(0.025,0.8,0.8,0.5,0.127,0.127,1,1,4,4,4,4,1,1)
#' }
#'
#'
#' @returns CP The consistency probability, a scalar.
#'
#'
conProb2MRCTs <- function(alpha,power1,power2,pi,rF1,rF2,d1,d2,sigmaTrt1,sigmaCtrl1,sigmaTrt2,sigmaCtrl2,randRatio1,randRatio2){
  N1 <- (randRatio1+1)*(sigmaTrt1^2/randRatio1+sigmaCtrl1^2)*(qnorm(1-alpha)+qnorm(power1))^2/(d1^2)
  N2 <- (randRatio2+1)*(sigmaTrt2^2/randRatio2+sigmaCtrl2^2)*(qnorm(1-alpha)+qnorm(power2))^2/(d2^2)
  sigmad1 <- sqrt((randRatio1+1)*(sigmaTrt1^2/randRatio1+sigmaCtrl1^2)/N1)
  sigmad2 <- sqrt((randRatio2+1)*(sigmaTrt2^2/randRatio2+sigmaCtrl2^2)/N2)
  w1 <- N1/(N1+N2)
  w2 <- N2/(N1+N2)
  integerand <- function(x){
    return(pnorm((1-pi)*(w1*sigmad1*x[1]+w2*sigmad2*x[2]+w1*d1+w2*d2)/sqrt((1/rF1-1)*w1^2*sigmad1^2+(1/rF2-1)*w2^2*sigmad2^2))*dnorm(x[1])*dnorm(x[2]))
  }
  CP <- cubature::adaptIntegrate(integerand,lower=c(-qnorm(power1),-qnorm(power2)),upper=c(Inf,Inf))
  CP <- CP$integral/(power1*power2)
  CP <- as.numeric(CP)
  return(CP)
}

