"_PACKAGE"

#' @title Regional consistency evaluation and sample size calculation for MRCTs
#'
#' @name ssMRCT-package
#' @description
#' The \pkg{ssMRCT} package offers a comprehensive suite of functions for designing and analyzing
#' Multi-Regional Clinical Trials (MRCTs), featuring specialized tools for:
#' \itemize{
#'   \item \strong{Consistency Probability Calculation}: Computes the probability that treatment effects
#'         across regions satisfy predefined consistency criteria;
#'   \item \strong{Optimal Regional Fraction Determination}: Identifies the ideal regional enrollment
#'         proportions to minimize regional sample size while maintaining consistency probability.
#' }
#' Support extends to both \emph{one MRCT} (one trial encompassing all regions) and \emph{two MRCTs}
#' (two pivotal, independent MRCTs) frameworks. Developed for regulatory decision
#' support and trial optimization in global drug development programs.
#'
#' @section Functions:
#' Key functions included in the package:
#' \describe{
#'   \item{\link{conProb}, \link{conProb2}}{
#'    Calculate the (conditional) consistency probability for one/two MRCT(s) 
#'    with corresponding overall sample size(s).
#'   }
#'   \item{\link{regFrac}, \link{regFrac2}}{
#'     Calculate the optimal regional fraction(s) given the (conditional) 
#'     consistency probability for one/two MRCT(s), with corresponding overall sample size(s).
#'   }
#' }
#' 
#' @references  
#' MHLW (2007). Basic Principles on Global Clinical Trials. https://www.pmda.go.jp/files/000153265.pdf
#' 
#' Kunhai Qing, Xinru Ren, Shuping Jiang, Ping Yang, Menggang Yu and Jin Xu (2025). Regional consistency evaluation and sample size calculation under two MRCTs. http://arxiv.org/abs/2411.15567
#' 
#'
#' @importFrom cubature adaptIntegrate
#' @importFrom stats dnorm pnorm qnorm sd integrate uniroot rbinom rnorm
NULL
