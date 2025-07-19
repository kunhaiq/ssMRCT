#' @name conProb2
#'
#' @title Consistency probability for two MRCTs
#'
#' @description Calculate the (conditional) consistency probability for two MRCTs.
#'
#' @param alpha The type I error.
#' @param power1 Power for MRCT 1.
#' @param power2 Power for MRCT 2. Defaults to \code{power1}.
#' @param pi The threshold ratio in Japan's criterion I (conditional version). Defaults to 0.5.
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
#' The extended (conditional) consistency probability, \eqn{\mathrm{Pr}\left(D_{k,\textrm{pool}}\ge \pi D_\textrm{pool}\ |\ T^{(1)}>z_{1-\alpha},T^{(2)}>z_{1-\alpha}\right)},
#' equals
#' \deqn{
#'  \begin{aligned}
#'  &\frac{1}{(1-\beta_{1})(1-\beta_{2})}\int_{-z_{1-\beta_{1}}}^{\infty}\int_{-z_{1-\beta_{2}}}^{\infty} \\
#'  &\quad
#'  \Phi\left(\frac{(1-\pi)\left\{w^{(1)}\sigma^{(1)}_{d}u+w^{(2)}\sigma^{(2)}_{d}v + w^{(1)}d^{(1)} + w^{(2)}d^{(2)}\right\}}{\sqrt{\left\{\left(f_{k}^{(1)}\right)^{-1}-1\right\} \left(w^{(1)}\sigma^{(1)}_{d}\right)^{2} + \left\{\left(f_{k}^{(2)}\right)^{-1}-1\right\} \left(w^{(2)}\sigma^{(2)}_{d}\right)^{2}}} \right)\phi(u)\phi(v)dudv.
#'  \end{aligned}
#' }
#' Assuming homogeneous variances, equal treatment effects
#' and equal randomization ratios across two studies, any of \code{d1},
#' \code{sigmaTrt1} and \code{randRatio1} is no required to provide. Then the
#' above equation reduces to
#' \deqn{
#'  \begin{aligned}
#'  &\frac{1}{(1-\beta_{1})(1-\beta_{2})}\int_{-(z_{1-\beta_{1}}+z_{1-\beta_{2}})/\sqrt{2}}^{\infty}\left[ \Phi\left(\frac{(1-\pi)\left\{\sqrt{2}u + (2z_{1-\alpha}+z_{1-\beta_{1}}+z_{1-\beta_{2}})\right\}}{\sqrt{\left(f_{k}^{(1)}\right)^{-1}+\left(f_{k}^{(2)}\right)^{-1}-2}}\right) \right. \\
#'  &\qquad\qquad\qquad\qquad \left.\left\{ \Phi\left(u+\sqrt{2}z_{1-\beta_{1}}\right)+\Phi\left(u+\sqrt{2}z_{1-\beta_{2}}\right)-1\right\} \right] \phi(u)du.
#'  \end{aligned}
#' }
#' Since there is no closed forms of above equations, \code{conProb2} utilizes
#' the \code{\link[cubature]{adaptIntegrate}} function or the \code{\link[stats]{integrate}}
#' function for numerical integration.
#'
#' @returns CP The consistency probability, a scalar.
#'
#' @examples
#' conProb2(0.025,0.8,0.8,0.5,0.127)
#'
#' conProb2(0.025,0.8,0.8,0.5,0.127,0.127,1,1,4,4,4,4,1,1)
#'
#' @export
#'
conProb2 <- function(alpha,power1,power2=NULL,pi=0.5,rF1,rF2=NULL,d1=NULL,d2=NULL,sigmaTrt1=NULL,sigmaCtrl1=NULL,sigmaTrt2=NULL,sigmaCtrl2=NULL,randRatio1=NULL,randRatio2=NULL){
  if (is.null(power2)) {
    power2 <- power1
  }
  if (is.null(rF2)) {
    rF2 <- rF1
  }
  if (is.null(sigmaTrt1) || is.null(d1) || is.null(randRatio1)) {
    cat("considering two homogeneous MRCTs \n")
    return(conProb2MRCTsH(alpha,power1,power2,pi,rF1,rF2))
  }
  if (is.null(sigmaCtrl1)) {
    sigmaCtrl1 <- sigmaTrt1
  }
  if (is.null(sigmaTrt2)) {
    sigmaTrt2 <- sigmaTrt1
  }
  if (is.null(sigmaCtrl2)) {
    sigmaCtrl2 <- sigmaTrt2
  }
  if (is.null(randRatio2)) {
    randRatio2 <- randRatio1
  }
  return(conProb2MRCTs(alpha,power1,power2,pi,rF1,rF2,d1,d2,sigmaTrt1,sigmaCtrl1,sigmaTrt2,sigmaCtrl2,randRatio1,randRatio2))
}

