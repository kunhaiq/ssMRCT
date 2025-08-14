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
#' @param sigmaTrt1 The standard deviation of response in the treatment group for MRCT 1.
#' @param sigmaCtrl1 The standard deviation of response in the control group for MRCT 1. Defaults to \code{sigmaTrt1}.
#' @param sigmaTrt2 The standard deviation of response in the treatment group for MRCT 2. Defaults to \code{sigmaTrt1}.
#' @param sigmaCtrl2 The standard deviation of response in the control group for MRCT 2. Defaults to \code{sigmaTrt2}.
#' @param randRatio1 The randomization ratio between the treatment group and control group for MRCT 1.
#' @param randRatio2 The randomization ratio between the treatment group and control group for MRCT 2. Defaults to \code{randRatio1}.
#'
#' @details
#' The extended consistency probability via the extended criterion Method I of MHLW (2007) (conditional version),
#' \deqn{\mathrm{Pr}\left(D_{k,\textrm{pool}}\ge \pi D_\textrm{pool}\ |\ T^{(1)}>z_{1-\alpha},T^{(2)}>z_{1-\alpha}\right),}
#' is approximately
#' \deqn{
#'  \begin{aligned}
#'  &\frac{1}{(1-\beta_{1})(1-\beta_{2})}\int_{-z_{1-\beta_{1}}}^{\infty}\int_{-z_{1-\beta_{2}}}^{\infty} \\
#'  &\quad
#'  \Phi\left(\frac{(1-\pi)\left\{w^{(1)}\sigma^{(1)}_{d}u+w^{(2)}\sigma^{(2)}_{d}v + w^{(1)}d^{(1)} + w^{(2)}d^{(2)}\right\}}{\sqrt{\left\{\left(f_{k}^{(1)}\right)^{-1}-1\right\} \left(w^{(1)}\sigma^{(1)}_{d}\right)^{2} + \left\{\left(f_{k}^{(2)}\right)^{-1}-1\right\} \left(w^{(2)}\sigma^{(2)}_{d}\right)^{2}}} \right)\phi(u)\phi(v)dudv,
#'  \end{aligned}
#' } 
#' where \eqn{w^{(s)}=N^{(s)}/(N^{(1)}+N^{(2)})} and 
#' \deqn{
#'  \sigma_d^{2}=\mathrm{var}(D)=\frac{(r+1)\left\{\sigma^{2(\textrm{t})}+r\sigma^{2(\textrm{c})}\right\}}{rN}.
#' }
#' 
#' Since there is no closed forms of above equations, \code{conProb2} utilizes
#' the \code{\link[cubature]{adaptIntegrate}} function for numerical integration.
#' 
#' Additionally, if we assume homogeneous variances, equal treatment effects
#' and equal randomization ratios across two studies, i.e., \eqn{\sigma^{2(\textrm{t},1)}=\sigma^{2(\textrm{t},2)}},
#' \eqn{\sigma^{2(\textrm{c},1)}=\sigma^{2(\textrm{c},2)}}, \eqn{d^{(1)}=d^{(2)}}
#' and \eqn{r^{(1)}=r^{(2)}}, then the above equation reduces to
#' \deqn{
#'  \begin{aligned}
#'  &\frac{1}{(1-\beta_{1})(1-\beta_{2})}\int_{-(z_{1-\beta_{1}}+z_{1-\beta_{2}})/\sqrt{2}}^{\infty}\left[ \Phi\left(\frac{(1-\pi)\left\{\sqrt{2}u + (2z_{1-\alpha}+z_{1-\beta_{1}}+z_{1-\beta_{2}})\right\}}{\sqrt{\left(f_{k}^{(1)}\right)^{-1}+\left(f_{k}^{(2)}\right)^{-1}-2}}\right) \right. \\
#'  &\qquad\qquad\qquad\qquad \left.\left\{ \Phi\left(u+\sqrt{2}z_{1-\beta_{1}}\right)+\Phi\left(u+\sqrt{2}z_{1-\beta_{2}}\right)-1\right\} \right] \phi(u)du.
#'  \end{aligned}
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
#' ### Remark 7
#' conProb2(alpha = 0.05, power1 = 0.8, power2 = 0.9, rF1 = 0.141, d1 = 1, sigmaTrt1 = 4)
#' conProb2(alpha = 0.05, power1 = 0.8, power2 = 0.9, rF1 = 0.100, rF2 = 0.238, d1 = 1, sigmaTrt1 = 4)
#' ### Remark 9
#' conProb2(alpha = 0.05, power1 = 0.8, rF1 = 0.154, d1 = 1, sigmaTrt1 = 4)
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
                     sigmaTrt1,
                     sigmaCtrl1 = sigmaTrt1,
                     sigmaTrt2 = sigmaTrt1,
                     sigmaCtrl2 = sigmaTrt2,
                     randRatio1 = 1,
                     randRatio2 = randRatio1){ 
  ### consistency probability
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
  ### sample size
  NTrt1 <-  ceiling((1/randRatio1*sigmaTrt1^2 + sigmaCtrl1^2)*(qnorm(1-alpha) + qnorm(power1))^2/d1^2)
  NCtrl1 <- ceiling(randRatio1*NTrt1)
  N1 <- NTrt1 + NCtrl1
  NTrt2 <-  ceiling((1/randRatio2*sigmaTrt2^2 + sigmaCtrl2^2)*(qnorm(1-alpha) + qnorm(power2))^2/d2^2)
  NCtrl2 <- ceiling(randRatio2*NTrt2)
  N2 <- NTrt2 + NCtrl2
  return(list(CP = CP, N1 = N1, N2 = N2)) 
}

