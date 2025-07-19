#' @name regFrac2
#'
#' @title Regional fraction for two MRCTs
#'
#' @description Calculate the optimal regional fractions given the (conditional) consistency probability for two MRCT.
#'
#' @param alpha The type I error.
#' @param power1 Power for MRCT 1.
#' @param power2 Power for MRCT 2. Defaults to \code{power1}.
#' @param pi The threshold ratio in Japan's criterion I (conditional version). Defaults to 0.5.
#' @param CP The consistency probability.
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
#' Given the (conditional) consistency probability, there is an optimal pair of regional
#' fractions \code{rF1} and \code{rF2} that minimized the combined region sample size,
#' i.e., \eqn{f_{k}^{(1)}N^{(1)} + f_{k}^{(2)}N^{(2)}}. Theoretically,
#' the ratio between such \code{rF1} and \code{rF2} is fixed. Hence,
#' \code{regFrac2} utilizes two core computational components:
#' \itemize{
#'   \item The \code{\link{conProb2}} function to compute the (conditional) consistency probability for two MRCTs
#'   \item The \code{\link[stats]{uniroot}} function from the \pkg{stats} package for numerical root-finding
#' }
#'
#' The solution is obtained by solving the following equation numerically:
#' \deqn{
#' \begin{array}{c}
#'   \code{conProb2(rF2*k,rF2)} - \code{CP} = 0 \\
#'   \text{where } \code{rF2} \in (0,\min(1/k,1)),
#' \end{array}
#' }
#' where \eqn{k} is the above fixed ratio between \code{rF1} and \code{rF2}.
#'
#' Assuming homogeneous variances, equal treatment effects
#' and equal randomization ratios across two studies, any of \code{d1},
#' \code{sigmaTrt1} and \code{randRatio1} is no required to provide.
#'
#' @return A list containing the following two components:
#' \describe{
#'   \item{\code{rF1}}{The optimal regional fraction for MRCT 1.}
#'   \item{\code{rF2}}{The optimal regional fraction for MRCT 2.}
#' }
#'
#' @examples
#'
#' regFrac2(0.025,0.8,0.8,0.5,0.8)
#'
#' regFrac2(0.025,0.8,0.8,0.5,0.8,1,1,4,4,4,4,1,1)
#'
#' @export
#'
regFrac2 <- function(alpha,power1,power2=NULL,pi=0.5,CP,d1=NULL,d2=NULL,sigmaTrt1=NULL,sigmaCtrl1=NULL,sigmaTrt2=NULL,sigmaCtrl2=NULL,randRatio1=NULL,randRatio2=NULL){
  if (is.null(power2)) {
    power2 <- power1
  }
  if (is.null(sigmaTrt1) || is.null(d1) || is.null(randRatio1)) {
    cat("considering two homogeneous MRCTs \n")
    return(regFrac2MRCTsH(alpha,power1,power2,pi,CP))
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
  return(regFrac2MRCTs(alpha,power1,power2,pi,CP,d1,d2,sigmaTrt1,sigmaCtrl1,sigmaTrt2,sigmaCtrl2,randRatio1,randRatio2))
}
