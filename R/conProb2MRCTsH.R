#' @name conProb2MRCTsH
#'
#' @title Consistency probability for two homogenous MRCTs
#'
#' @description Calculate the (conditional) consistency probability for two homogeneous MRCTs
#'
#' @param alpha The type I error.
#' @param power1 Power for MRCT 1.
#' @param power2 Power for MRCT 2.
#' @param pi The threshold ratio in Japan's criterion I (conditional version).
#' @param rF1 The regional fraction for MRCT 1.
#' @param rF2 The regional fraction for MRCT 2.
#'
#' @details
#' The extended (conditional) consistency probability, \eqn{\mathrm{Pr}\left(D_{k,\textrm{pool}}\ge \pi D_\textrm{pool}\ |\ T^{(1)}>z_{1-\alpha},T^{(2)}>z_{1-\alpha}\right)},
#' equals
#' \deqn{
#'  \begin{aligned}
#'  &\frac{1}{(1-\beta_{1})(1-\beta_{2})}\int_{-(z_{1-\beta_{1}}+z_{1-\beta_{2}})/\sqrt{2}}^{\infty}\left[ \Phi\left(\frac{(1-\pi)\left\{\sqrt{2}u + (2z_{1-\alpha}+z_{1-\beta_{1}}+z_{1-\beta_{2}})\right\}}{\sqrt{\left(f_{k}^{(1)}\right)^{-1}+\left(f_{k}^{(2)}\right)^{-1}-2}}\right) \right. \\
#'  &\qquad\qquad\qquad\qquad \left.\left\{ \Phi\left(u+\sqrt{2}z_{1-\beta_{1}}\right)+\Phi\left(u+\sqrt{2}z_{1-\beta_{2}}\right)-1\right\} \right] \phi(u)du,
#'  \end{aligned}
#' }
#' when assuming homogeneous variances, equal treatment effects
#' and equal randomization ratios across two studies,
#' Since there is no closed form of above equation, \code{conProb2MRCTsH} utilizes
#' the \code{\link[stats]{integrate}} function for numerical integration.
#'
#' @returns CP The consistency probability, a scalar.
#'
#' @examples
#' \dontrun{
#' conProb2MRCTsH(0.025,0.8,0.8,0.5,0.127,0.127)
#' }
#'
conProb2MRCTsH <- function(alpha,power1,power2,pi,rF1,rF2){
  integerand <- function(x){
    return(pnorm((1-pi)*(sqrt(2)*x+2*qnorm(1-alpha)+qnorm(power1)+qnorm(power2))/sqrt(1/rF1+1/rF2-2))*(pnorm(x+sqrt(2)*qnorm(power1))+pnorm(x+sqrt(2)*qnorm(power2))-1)*dnorm(x))
  }
  CP <- stats::integrate(integerand,lower=-(qnorm(power1)+qnorm(power2))/sqrt(2),upper=Inf)
  CP <- CP$value/(power1*power2)
  CP <- as.numeric(CP)
  return(CP)
}



