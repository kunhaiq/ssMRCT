#' @name regFrac2
#'
#' @title Regional fraction for two MRCTs
#'
#' @description Calculate the optimal regional fractions given the (conditional) consistency probability for two MRCTs.
#'
#' @param alpha The type I error.
#' @param power1 Power for MRCT 1.
#' @param power2 Power for MRCT 2. Defaults to \code{power1}.
#' @param pi The threshold ratio in Japan's criterion I (conditional version). Defaults to 0.5.
#' @param CP The consistency probability. Defaults to 80%.
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
#' The overall sample size is obtain by \code{\link{conProb2}} simultaneously. 
#' 
#' @returns A list containing the following two components:
#' \describe{
#'   \item{\code{rF1}}{The regional fraction for MRCT 1, a scalar.}
#'   \item{\code{rF2}}{The regional fraction for MRCT 2, a scalar.}
#'   \item{\code{N1}}{The overall sample size for MRCT 1, an integer.}
#'   \item{\code{N2}}{The overall sample size for MRCT 2, an integer.}
#' }
#' 
#' @examples 
#' 
#' regFrac2(alpha = 0.025, power1 = 0.8, d1 = 1, sigmaTrt1 = 4)
#'
#' @export
#'
regFrac2 <- function(alpha,
                     power1,
                     power2 = power1,
                     pi = 0.5,
                     CP = 0.8,
                     d1,
                     d2 = d1,
                     sigmaTrt1,
                     sigmaCtrl1 = sigmaTrt1,
                     sigmaTrt2 = sigmaTrt1,
                     sigmaCtrl2 = sigmaTrt2,
                     randRatio1 = 1,
                     randRatio2 = randRatio1){ 
  ### regional fraction
  N1 <- (randRatio1+1)*(sigmaTrt1^2/randRatio1+sigmaCtrl1^2)*(qnorm(1-alpha)+qnorm(power1))^2/(d1^2)
  N2 <- (randRatio2+1)*(sigmaTrt2^2/randRatio2+sigmaCtrl2^2)*(qnorm(1-alpha)+qnorm(power2))^2/(d2^2)
  sigmad1 <- sqrt((randRatio1+1)*(sigmaTrt1^2/randRatio1+sigmaCtrl1^2)/N1)
  sigmad2 <- sqrt((randRatio2+1)*(sigmaTrt2^2/randRatio2+sigmaCtrl2^2)/N2)
  k <- sigmad1*sqrt(N1)/sigmad2/sqrt(N2)
  f <- function(x){
    return(conProb2(alpha, power1, power2, pi, rF1 = k*x, rF2 = x, d1, d2, 
                    sigmaTrt1, sigmaCtrl1, sigmaTrt2, sigmaCtrl2, randRatio1, 
                    randRatio2)$CP - CP)
  }
  rF <- uniroot(f,interval=c(0,min(1/k,1)))$root
  rF <- list(rF1=k*rF,
             rF2=rF)
  ### sample size 
  N <- conProb2(alpha, power1, power2, pi, rF1 = rF$rF1, rF2 = rF$rF2, d1, d2, 
                 sigmaTrt1, sigmaCtrl1, sigmaTrt2, sigmaCtrl2, randRatio1, 
                 randRatio2)
  N1 <- N$N1
  N2 <- N$N2 
  return(list(rF1 = rF$rF1,
              rF2 = rF$rF2,
              N1 = N1,
              N2 = N2))
}
