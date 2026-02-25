#' @name regFrac2
#'
#' @title Regional fraction for two MRCTs via the extended criterion Method I of MHLW (2007) (conditional version) under fixed effects model
#'
#' @description Calculate the optimal regional fractions given the consistency probability for two MRCTs via the extended criterion Method I of MHLW (2007)  (conditional version) under fixed effects model
#'
#' @param alpha The Type I error.
#' @param power Power for MRCT 1 and 2. 
#' @param pi The threshold ratio in the extended criterion Method I of MHLW (2007) (conditional version). Defaults to 0.5.
#' @param CP The consistency probability. Defaults to 80%.
#'
#' @details
#' If we assume equal treatment effects, common regional fractions, common powers and common randomization ratios across 
#' two studies, i.e., \eqn{\sigma^{2(\textrm{t},1)}=\sigma^{2(\textrm{t},2)}, \sigma^{2(\textrm{c},1)}=\sigma^{2(\textrm{c},2)},
#' d^{(1)}=d^{(2)}, f_k^{(1)}=f_k^{(2)}=f_k, \beta_1=\beta_2} and \eqn{r^{(1)}=r^{(2)}}, then the CP reduces to
#' \deqn{
#'  \frac{1}{(1-\beta)^2} \int_{-z_{1-\beta}}^{\infty} \int_{-z_{1-\beta}}^{\infty} \Phi\left( \frac{(1-\pi) (2z_{1-\alpha} + 2z_{1-\beta} + u + v)}{\sqrt{2(f_k^{-1} - 1)}} \right) \phi(u)\phi(v) \mathrm{d}u \mathrm{d}v.
#' } 
#' Given the consistency probability, there is a minimal regional fraction \code{rF}. 
#'
#' The overall sample size is obtain by \code{\link{conProb2}} simultaneously. 
#' 
#' @returns A list containing the following one component:
#' \describe{
#'   \item{\code{rF}}{The regional fraction for MRCT 1 and 2, a scalar.}
#' }
#' 
#' @examples 
#' regFrac2(alpha = 0.05, power = 0.8)  
#' 
#' @export
#'
regFrac2 <-  function(alpha, 
                      power, 
                      pi = 0.5, 
                      CP = 0.8) { 
  za <- qnorm(1 - alpha)
  zb <- qnorm(power) 
  calculate_cp <- function(fk) { 
    if (fk <= 0 || fk >= 1) return(NA) 
    denom_val <- sqrt(2 * (1/fk - 1)) 
    integrand_func <- function(v) {
      sapply(v, function(vi) {
        
        inner_func <- function(u) { 
          numerator <- (1 - pi) * (2 * za + 2 * zb + u + vi) 
          pnorm(numerator / denom_val) * dnorm(u)
        } 
        res_u <- integrate(inner_func, lower = -zb, upper = Inf)$value 
        return(res_u * dnorm(vi))
      })
    } 
    integral_res <- integrate(integrand_func, lower = -zb, upper = Inf)$value 
    cp_val <- integral_res / ((power)^2)
    return(cp_val)
  }
   
  obj_func <- function(fk) {
    calculate_cp(fk) - CP
  } 
  root_res <- uniroot(obj_func, interval = c(0.001, 0.999), tol = 1e-6)
  return(root_res$root)
}
