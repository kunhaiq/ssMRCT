#' @name regFrac2MRCTsH
#'
#' @title Regional fraction for two homogeneous MRCTs
#'
#' @description Calculate the optimal regional fractions given the (conditional) consistency probability for two homogeneous MRCT.
#'
#' @param alpha The type I error.
#' @param power1 Power for MRCT 1.
#' @param power2 Power for MRCT 2.
#' @param pi The threshold ratio in Japan's criterion I (conditional version).
#' @param CP The consistency probability.
#'
#'
#' @details
#' Given the (conditional) consistency probability, there is an optimal pair of regional
#' fractions \code{rF1} and \code{rF2} that minimized the combined region sample size,
#' i.e., \eqn{f_{k}^{(1)}N^{(1)} + f_{k}^{(2)}N^{(2)}}. Theoretically,
#' the ratio between such \code{rF1} and \code{rF2} is fixed. Hence,
#' \code{regFrac2MRCTs} utilizes two core computational components:
#' \itemize{
#'   \item The \code{\link{conProb2MRCTs}} function to compute the (conditional) consistency probability for two MRCTs
#'   \item The \code{\link[stats]{uniroot}} function from the \pkg{stats} package for numerical root-finding
#' }
#' The solution is obtained by solving the following equation numerically:
#' \deqn{
#' \begin{array}{c}
#'   \code{conProb2MRCTs(rF2*k,rF2)} - \code{CP} = 0 \\
#'   \text{where } \code{rF2} \in (0,\min(1/k,1)),
#' \end{array}
#' }
#' where \eqn{k} is the above fixed ratio between \code{rF1} and \code{rF2}.
#'
#' @return A list containing the following two components:
#' \describe{
#'   \item{\code{rF1}}{The optimal regional fraction for MRCT 1.}
#'   \item{\code{rF2}}{The optimal regional fraction for MRCT 2.}
#' }
#'
#' @examples
#' \dontrun{
#' regFrac2MRCTsH(0.025,0.8,0.8,0.5,0.8)
#' }
#'
regFrac2MRCTsH <- function(alpha,power1,power2,pi,CP){
  k <- 1
  f <- function(x){
    return(conProb2MRCTsH(alpha,power1,power2,pi,rF1=k*x,rF2=x)-CP)
  }
  rF <- uniroot(f,interval=c(0,min(1/k,1)))$root
  rF <- list(rF1=k*rF,
             rF2=rF)
  return(rF)
}
