#' @name regFrac
#'
#' @title Regional fraction for one MRCT
#'
#' @description Calculate the optimal regional fraction given the (conditional) consistency probability for one MRCT.
#'
#' @param alpha The type I error.
#' @param power Power.
#' @param pi The threshold ratio in Japan's criterion I (conditional version). Defaults to 0.5.
#' @param CP the (conditional) consistency probability.
#'
#' @details
#' Given the (conditional) consistency probability, there is a minimum regional
#' fraction \code{rF}. To calculate the minimum \code{rF}, \code{regFrac} utilizes two core computational components:
#' \itemize{
#'   \item The \code{\link{conProb}} function to compute the (conditional) consistency probability
#'   \item The \code{\link[stats]{uniroot}} function from the \pkg{stats} package for numerical root-finding
#' }
#' The solution is obtained by solving the following equation numerically:
#' \deqn{
#' \begin{array}{c}
#'   \code{conProb(rF)} - \code{CP} = 0 \\
#'   \text{where } \code{rF} \in (0,1).
#' \end{array}
#' }
#' @returns rF The regional fraction, a scalar.
#'
#' @examples
#'
#' regFrac(alpha=0.025,power=0.8,pi=0.5,CP=0.8)
#'
#' @export
#'
regFrac <- function(alpha,power,pi=0.5,CP){
  f <- function(x){
    return(conProb(alpha,power,pi,rF = x)-CP)
  }
  rF <- stats::uniroot(f,interval = c(0,1))$root
  rF <- as.numeric(rF)
  return(rF)
}
