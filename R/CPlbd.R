#' @name CPlbd
#'
#' @title Lower bound of the consistency probability for one MRCT via Method I of MHLW (2007)(conditional version) under random effects model.
#'
#' @description Calculate the lower bound of consistency probability for one MRCT via Method I of MHLW (2007)(conditional version) under random effects model for normal endpoint.
#'
#' @param alpha The one-sided Type I error.
#' @param beta The Type II error.
#' @param tau_delta_r The ratio of the standard deviation vs mean of the Gaussian prior for the regional treatment effect.
#' @param delta The mean of the Gaussian prior for the regional treatment effect.
#' @param sigma0r The standard deviation(s) of response in the treatment group for region(s) of interest. A scalar or vector of length equals to the total number of region(s) of interest.
#' @param sigma1r The standard deviation(s) of response in the control group for region(s) of interest. A scalar or vector of length equals to \code{sigma0r}.
#' @param pi The threshold ratio in Method I of MHLW (2007)(conditional version). Defaults to 0.5.
#' @param randRatio The randomization ratio between the treatment group and control group. Defaults to 1.
#' @param fr The regional fraction(s) for the region(s) of interest. A scalar or vector of equal length to \code{sigma0r} and \code{sigma1r}.
#'
#'
#' @details
#' The lower bound of the consistency probability via Method I of MHLW (2007)(conditional version) under random effects model,
#' \eqn{ \mathrm{Pr} \left( \widetilde{D}_{r} \ge \pi \widetilde{D} | T > z_{1-\alpha} \right) },
#' is approximately
#' \deqn{
#'  \frac{1}{1-\beta}\int_{u>z_{\beta}}\Phi\left(\frac{2\delta^2(1-\pi)(u+z_{1-\alpha}+z_{1-\beta})}{\tau^2(z_{1-\alpha}+z_{1-\beta})^2} \right) d\Phi(u).
#' }
#' provided \eqn{\frac{\tau}{\delta}<\frac{\sqrt{2}}{z_{1-\alpha}+z_{1-\beta}}}. Otherwise,
#' CP will also be greater than
#' \deqn{
#' \frac{1}{1-\beta}\int_{u>z_{\beta}}\Phi\left(\frac{(1-\pi) \left(u+z_{1-\alpha}+z_{1-\beta}\right) }{\sqrt{\frac{\tau^2}{\delta^2}(z_{1-\alpha}+z_{1-\beta})^2-1}}\right) d\Phi(u).
#' }.
#' Since there is no closed form of above equation, \code{CPlbd} utilizes the \code{\link[stats]{integrate}} function for numerical integration.
#'
#' The sample size(s) for the region(s) of interest at the lower bound are calculated based on solving the following equation for \eqn{n_0^{(r)}},
#' \deqn{
#' \frac{h_r}{h_r+1}=\frac{\tau^2}{2\delta^2}(z_{1-\alpha}+z_{1-\beta})^2, \quad n_1^{(r)}=\ell n_0^{(r)}
#' }
#' where \eqn{h_r=\tau^2/\sigma^{2(r)}, \sigma^{2(r)}=\Omega_r/(n_0f_r), \Omega_r=\ell^{-1}\sigma_1^{2(r)}+\sigma_1^{2(r)}}.
#' Then \eqn{n^{(r)} = n_0^{(r)}+n_1^{(r)}}, where both of
#' \eqn{n_0^{(r)}} and \eqn{n_1^{(r)}} should be integers and hence \eqn{n^{(r)}}.
#'
#' @returns A list containing the following components:
#' \describe{
#'   \item{\code{CPlbd}}{The lower bound of consistency probability, a scalar.}
#'   \item{\code{n0r}}{The regional sample size of the control group. (n.a. when the condition for positive sample size is not met)}
#'   \item{\code{n1r}}{The regional sample size of the treatment group. (n.a. when the condition for positive sample size is not met)}
#'   \item{\code{nr}}{The regional sample size. (n.a. when the condition for positive sample size is not met)}
#'
#'
#' }
#'
#' @examples
#' CPlbd(alpha=0.05, beta=0.2, tau_delta_r=0.4, pi=0.5, delta=0.25,
#'       sigma0r=1, sigma1r=1, randRatio=1, fr=c(0.1))
#' CPlbd(alpha=0.05, beta=0.2, tau_delta_r=0.4, pi=0.5, delta=0.25,
#'       sigma0r=c(1,1), sigma1r=c(1,1), randRatio=1, fr=c(0.1,0.5))
#'
#' @export
#'
CPlbd<-function(alpha, beta, tau_delta_r, delta, sigma0r,sigma1r,fr,pi=0.5,randRatio=1){  #tau_delta_r is the ratio of tau divided by delta
  threshold<-sqrt(2)/(qnorm(1-alpha)+qnorm(1-beta))  # the condition to determine the right cp lower bound eq to use

  ### condition met
  if(tau_delta_r<threshold){   # the "when" condition of corollary 1
    integrand1<-function(x){
      return(pnorm(2*tau_delta_r^(-2)*(1-pi)*(x+qnorm(1-beta)+qnorm(1-alpha))/(qnorm(1-beta)+qnorm(1-alpha))^2)*dnorm(x))  # Phi(.) to be integrated in cp lower bound expression in eq(13)
    }
    out1<-stats::integrate(integrand1, lower=qnorm(beta), upper=Inf)  # perform integration
    cpmin<-out1$value/(1-beta)   # cp lower bound

    tau<-tau_delta_r*delta
    a<-tau_delta_r^2*(qnorm(1-alpha)+qnorm(1-beta))^2/2    # RHS on the '=' of the when condition in eq (13)
    Omega<-sigma0r^2+randRatio^(-1)*sigma1r^2
    n0r<-ceiling(Omega/((1/a-1)*fr*tau^2))    # expression of n0 by solving the "=" of the when condition
    n1r<-ceiling(n0r*randRatio);
    nr<-n0r+n1r
  }


  ### condition not met
  else {
    integrand2<-function(x){
      return(pnorm((1-pi)*(x+qnorm(1-beta)+qnorm(1-alpha))/sqrt((tau_delta_r*(qnorm(1-alpha)+qnorm(1-beta)))^2-1))*dnorm(x))  # Phi(.) to be integrated in cp lower bound expression in eq(14)
    }
    out2<-stats::integrate(integrand2, lower=qnorm(beta), upper=Inf)  # perform integration
    cpmin<-out2$value/(1-beta)   # cp lower bound per (14)
    n0r<-"n.a.";n1r<-"n.a."; nr<-"n.a."
  }

  return(list("CPlbd"=round(cpmin,3), "n0r"=n0r,"n1r"=n1r,"nr"=nr))
}


