#' @name conProb
#'
#' @title Consistency probability for one MRCT
#'
#' @description Calculate the (conditional) consistency probability for one MRCT.
#'
#' @param alpha The Type I error.
#' @param power Power.
#' @param pi The threshold ratio in Japan's criterion I (conditional version). Defaults to 0.5.
#' @param rF The regional fraction.
#'
#' @details
#' The (conditional) consistency probability, \eqn{\mathrm{Pr}\left(D_{k}\ge \pi D\ |\ T>z_{1-\alpha}\right)},
#' equals
#' \deqn{
#'  \frac{1}{1-\beta}\int_{-z_{1-\beta}}^{\infty} \Phi\left(\frac{(1-\pi)(u + z_{1-\alpha}+z_{1-\beta})}{\sqrt{f_{k}^{-1}-1}}\right)\phi(u)du.
#' }
#' Since there is no closed form of above equation, \code{conProb} utilizes the \code{\link[stats]{integrate}} function for numerical integration.
#'
#' @returns CP The consistency probability, a scalar.
#'
#' @examples
#' conProb(alpha=0.025,power=0.8,pi=0.5,rF=0.23)
#'
#' @export
#'
conProb <- function(alpha,power,pi=0.5,rF){
  integerand <- function(x){
    return(pnorm((1-pi)*(x+qnorm(1-alpha)+qnorm(power))/sqrt(1/rF-1))*dnorm(x))
  }
  CP <- stats::integrate(integerand,lower=-qnorm(power),upper=Inf)
  CP <- CP$value/power
  CP <- as.numeric(CP)
  return(CP)
}
