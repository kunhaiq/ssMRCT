"_PACKAGE"

#' @title Regional consistency evaluation and sample size calculation for MRCTs
#'
#' @name ssMRCT-package
#' @description
#' The \pkg{ssMRCT} package offers a comprehensive suite of functions for designing and analyzing
#' Multi-Regional Clinical Trials (MRCTs), featuring specialized tools for:
#' \itemize{
#'   \item \strong{Consistency Probability Calculation}: Computes the probability that treatment effects
#'         across regions satisfy predefined consistency criteria
#'          under fixed effects model and random effects model;
#'   \item \strong{Optimal Regional Fraction Determination}: Identifies the ideal regional enrollment
#'         proportions to minimize regional sample size while maintaining consistency probability
#'         under fixed effects model.
#' }
#' Support extends two MRCTs (two pivotal, independent MRCTs) under fixed effects model. 
#' Developed for regulatory decision support and trial optimization in global drug development programs.
#' In this package, all sample sizes are rounded up to be integers. 
#' @section Functions:
#' Key functions included in the package:
#' \describe{
#'   \item{\link{conProb}, \link{conProb2}}{
#'    Calculate the (conditional) consistency probability under criterion Method I of MHLW (2007)
#'    for one/two MRCT(s) with corresponding overall sample size(s) under fixed effects model.
#'   }
#'   \item{\link{regFrac}, \link{regFrac2}}{
#'     Calculate the optimal regional fraction(s) given the (conditional) 
#'     consistency probability under criterion Method I of MHLW (2007)
#'     for one/two MRCT(s), with corresponding overall sample size(s) under fixed effects model.
#'   }
#'   \item{\link{ssCPnorm}, \link{ssCPbinary}}{
#'     Calculate the (conditional) consistency probability under criterion Method I of MHLW (2007)
#'     for MRCT with corresponding overall sample size under random effects model
#'     with normal and binary response, respectively.
#'   }
#'   \item{\link{ssCPsurvPH}, \link{ssCPsurvNPH}}{
#'     Calculate the (conditional) consistency probability under criterion Method I of MHLW (2007)
#'     for MRCT with corresponding overall sample size under random effects model
#'     with survival endpoints (proportional or non-proportional hazard).
#'   }
#' }
#' 
#' @references  
#' MHLW (2007). Basic Principles on Global Clinical Trials. https://www.pmda.go.jp/files/000153265.pdf
#' 
#' Kunhai Qing, Xinru Ren, Shuping Jiang, Ping Yang, Menggang Yu and Jin Xu (2025). Regional consistency evaluation and sample size calculation under two MRCTs. http://arxiv.org/abs/2411.15567
#' 
#' Xinru Ren, Jin Xu (2025). Consistency assessment and regional sample size calculation for MRCTs under random effects model.  https://arxiv.org/abs/2508.09443
#'
#' @importFrom cubature adaptIntegrate
#' @importFrom stats dnorm pnorm qnorm sd integrate uniroot rbinom rnorm
NULL
