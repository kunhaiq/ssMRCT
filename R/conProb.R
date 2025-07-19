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
#' @param d The true mean of difference of response.
#' @param sigmaTrt The standard deviation of response in the treatment group.
#' @param sigmaCtrl The standard deviation of response in the control group. Defaults to \code{sigmaTrt}.
#' @param randRatio The randomization ratio between the treatment group and control group. Defaults to 1.
#'
#' @details
#' The (conditional) consistency probability, \eqn{\mathrm{Pr}\left(D_{k}\ge \pi D\ |\ T>z_{1-\alpha}\right)},
#' equals
#' \deqn{
#'  \frac{1}{1-\beta}\int_{-z_{1-\beta}}^{\infty} \Phi\left(\frac{(1-\pi)(u + z_{1-\alpha}+z_{1-\beta})}{\sqrt{f_{k}^{-1}-1}}\right)\phi(u)du.
#' }
#' Since there is no closed form of above equation, \code{conProb} utilizes the \code{\link[stats]{integrate}} function for numerical integration.
#'   
#' As the first equation includes none of \code{d}, \code{sigmaTrt}, 
#' \code{sigmaCtrl} and \code{randRatio}, if only the consistency probability
#' is considered, then the values of \code{d} and \code{sigmaTrt} could be arbitrary.
#' 
#' The overall sample size is calculated based on the below equation,
#' \deqn{
#' N^{(\textrm{c})}=\frac{\left\{r^{-1}\sigma^{2(\textrm{t})}+\sigma^{2(\textrm{c})}\right\}(z_{1-\alpha}+z_{1-\beta})^{2}}{d^{2}}, \quad N^{(\textrm{t})}=rN^{(\textrm{c})}.
#' }
#' Then \eqn{N = N^{(\textrm{t})} + N^{(\textrm{c})}}. Additionally, both of 
#' \eqn{N^{(\textrm{t})}} and \eqn{N^{(\textrm{c})}} should be integers and hence \eqn{N}.
#' 
#' @returns A list containing the following two components:
#' \describe{
#'   \item{\code{CP}}{The consistency probability, a scalar.}
#'   \item{\code{N}}{The overall sample size.}
#' } 
#' 
#' @examples
#' 
#' conProb(alpha = 0.025, power = 0.8, rF = 0.23, d = 1, sigmaTrt = 4)
#'
#' @export
#'
conProb <- function(alpha, 
                    power, 
                    pi = 0.5, 
                    rF, 
                    d, 
                    sigmaTrt, 
                    sigmaCtrl = sigmaTrt, 
                    randRatio = 1){
  ### consistency probability
  integerand <- function(x){
    return(pnorm((1-pi)*(x+qnorm(1-alpha)+qnorm(power))/sqrt(1/rF-1))*dnorm(x))
  }
  CP <- stats::integrate(integerand,lower=-qnorm(power),upper=Inf)
  CP <- CP$value/power
  CP <- as.numeric(CP)  
  ### sample size
  NTrt <-  ceiling((1/randRatio*sigmaTrt^2 + sigmaCtrl^2)*(qnorm(1-alpha) + qnorm(power))^2/d^2)
  NCtrl <- ceiling(randRatio*NTrt)
  N <- NTrt + NCtrl
  return(list(CP = CP, N = N))
}
