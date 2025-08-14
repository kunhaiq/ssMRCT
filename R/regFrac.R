#' @name regFrac
#'
#' @title Regional fraction for MRCT via criterion Method I of MHLW (2007) (conditional version) under fixed effects model
#'
#' @description Calculate the minimal regional fraction given the consistency probability for MRCT via criterion Method I of MHLW (2007) (conditional version) under fixed effects model.
#'
#' @param alpha The Type I error.
#' @param power Power.
#' @param pi The threshold ratio in criterion Method I of MHLW (2007) (conditional version). Defaults to 0.5.
#' @param CP the consistency probability. Defaults to 80%.
#' @param d The true mean of difference of response.
#' @param sigmaTrt The standard deviation of response in the treatment group.
#' @param sigmaCtrl The standard deviation of response in the control group. Defaults to \code{sigmaTrt}.
#' @param randRatio The randomization ratio between the treatment group and control group. Defaults to 1.
#'
#' @details
#' Given the consistency probability, there is a minimal regional
#' fraction \code{rF}. To calculate the minimal \code{rF}, \code{regFrac} utilizes two core computational components:
#' \itemize{
#'   \item The \code{\link{conProb}} function to compute the (conditional) consistency probability
#'   \item The \code{\link[stats]{uniroot}} function from the \pkg{stats} package for numerical root-finding
#' }
#' 
#' The solution is obtained by solving the following equation numerically:
#' \deqn{
#' \begin{array}{c}
#'   \code{conProb(rF)} - \code{CP} = 0 \\
#'   \text{where } \code{rF} \in (0,1).
#' \end{array}
#' }
#' 
#' The overall sample size is obtain by \code{\link{conProb}} simultaneously.
#' 
#' @returns A list containing the following two components:
#' \describe{
#'   \item{\code{rF}}{The regional fraction, a scalar.}
#'   \item{\code{N}}{The overall sample size.}
#' } 
#' 
#'
#' @examples
#' conProb(alpha = 0.05, power = 0.8, rF = 0.271, d = 1, sigmaTrt = 4)
#' 
#' regFrac(alpha = 0.05, power = 0.8, d = 1, sigmaTrt = 4)
#'   
#' @export
#'
regFrac <- function(alpha, 
                    power, 
                    pi = 0.5, 
                    CP = 0.8, 
                    d, 
                    sigmaTrt, 
                    sigmaCtrl = sigmaTrt, 
                    randRatio = 1){
  ### regional fraction
  f <- function(x){
    return(conProb(alpha, power, pi, rF = x, d, sigmaTrt, sigmaCtrl, randRatio)$CP-CP)
  }
  rF <- stats::uniroot(f,interval = c(0,1))$root
  rF <- as.numeric(rF)
  ### sample size
  N <- conProb(alpha, power, pi, rF, d, sigmaTrt, sigmaCtrl, randRatio)$N 
  return(list(rF = rF, N = N))
}
