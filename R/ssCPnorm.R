#' @name ssCPnorm
#'
#' @title Overall sample size and consistency probability for one MRCT via Method I of MHLW (2007)(conditional version) under random effects model.
#'
#' @description Calculate the overall sample size and consistency probability for one MRCT via Method I of MHLW (2007)(conditional version) under random effects model for the normal endpoint.
#'
#' @param alpha Type I error for one-sided test.
#' @param beta Type II error.
#' @param delta mean of the Gaussian prior for the regional treatment effect.
#' @param tau standard deviation of the Gaussian prior for the regional treatment effect.
#' @param sigma1r regional standard deviations of the response in the treatment group, a \code{R x 1} vector with each component representing the standard deviation in region \eqn{r} .
#' @param sigma0r regional standard deviations of the response in the control group. a \code{R x 1} vector with each component representing the standard deviation in region \eqn{r}.
#' @param pi threshold ratio in Method I of MHLW (2007)(conditional version). Defaults to 0.5.
#' @param randRatio randomization ratio between the treatment group and control group. Defaults to 1.
#' @param target The index for the region of interest. an integer of value referring to the region with the corresponding position in \code{f}.
#' @param f regional fractions, a \code{R x 1} vector with each component representing the fraction in region \eqn{r} and the sum must equal 1.
#'
#'
#' @details
#' The overall sample size is calculated based on the following equation,
#' \deqn{
#' \sum_r\left(\tau^2+\frac{\Omega_r}{n_0f_r}\right)^{-1} =\frac{(z_{1-\alpha}+z_{1-\beta})^2}{\delta^2} \quad n_1=\ell n_0.
#' }
#' Since there is no closed form of above equation, \code{ssCPnorm} utilizes the \code{\link[stats]{uniroot}} function for numerical solution for \eqn{n_0}.
#' Then \eqn{n = n_0+n_1}.
#'
#' The consistency probability via Method I of MHLW (2007)(conditional version) under random effects model,
#' \eqn{\mathrm{Pr}\left(\widetilde{D}_{r}\ge \pi \widetilde{D} | T>z_{1-\alpha}\right)},
#' is approximately
#' \deqn{
#'  \frac{1}{1-\beta}\int_{u>z_{\beta}}\Phi\left(\frac{(1-\pi) \left(u+z_{1-\alpha}+z_{1-\beta}\right) }{\sqrt{\rho_r^{-1}-1}}\right) d\Phi(u)
#' }
#' where \eqn{\rho_r^{-1}-1=\frac{h_r}{h_r+1}\sum_{j \ne r}^R\frac{h_j}{h_j+1},
#' w=\frac{1}{\tau^2}\sum_{j=1}^R\frac{h_j}{h_j+1}, h_r=\tau^2/\sigma^{2(r)}, \sigma^{2(r)}=\Omega_r/(n_0f_r)}.
#' Since there is no closed form of above equation, \code{ssCPnorm} utilizes the \code{\link[stats]{integrate}} function for numerical integration.
#'
#'For the normal endpoint, \eqn{\Omega_r=\ell^{-1}\sigma_1^{2(r)}+\sigma_1^{2(r)}}
#'
#' @returns A list containing the following components:
#' \describe{
#'   \item{\code{n0}}{The overall sample size of the control group}
#'   \item{\code{n1}}{The overall sample size of the treatment group}
#'   \item{\code{n}}{The overall sample size.}
#'   \item{\code{CP}}{The consistency probability, a scalar.}
#'   \item{\code{"error: condition not met. negative n0 returned"}}{}
#'   \item{\code{"error: regional fractions do not sum up to 1"}}{When vector f does not sum up to 1}
#' }
#'
#' @examples
#'
#' ssCPnorm(alpha=0.025, beta=0.1, delta=0.25, tau=0.1,
#'          sigma1r=rep(1,3),sigma0r=rep(1,3),
#'          randRatio=1, pi=0.5, target=1, f=c(1/3,1/3,1/3))
#'
#' @export
#'
ssCPnorm<-function(alpha, beta, delta, tau, sigma1r, sigma0r, pi=0.5, randRatio=1, target, f){

  R<-length(f);
  tau_delta_r<-tau/delta
  condition<-sqrt(R)/(qnorm(1-alpha)+qnorm(1-beta))  ## upper value of tau_delta_r (for homogeneous Omegars and f)

  if(abs(sum(f)-1) > 10^(-10)){   ## instead of if(sum(f)!=1) -- to Use a tolerance-based comparison for floating-point numbers, which introduce tiny rounding errors (e.g. (1-0.6)==(3*0.28) returns FALSE)
    return(list("error: regional fractions do not sum up to 1"))
    }
  else if (tau_delta_r<condition){ #### add a check to ensure positive n0 to fit into the context
    Omegars<-randRatio^(-1)*sigma1r^2+sigma0r^2

    #### calculate overall sample size n0
    fp<-function(n0){
      rhs <- (qnorm(1-alpha)+qnorm(1-beta))^2/delta^2
      sum((tau^2+Omegars/(n0*f))^(-1))-rhs
    }
    n0<-ceiling(stats::uniroot(fp,c(1,100000), tol = 0.0001)$root) 	# solve eq.
    n1<-ceiling(n0*randRatio)
    #print(n0)

    #### step 2 calculating CP using formula
    integrand<-function(x){
      b<-1/(1+Omegars[target]/(tau^2*n0*f[target]))   # replace f with f[1] to accommodate different fs
      c<-sum(1/(1+Omegars[-target]/(tau^2*n0*f[-target]))) # same as above f[-1]
      deno=sqrt(b*c)
      pnorm((1-pi)*(x+qnorm(1-beta)+qnorm(1-alpha))/deno)*dnorm(x)
    }
    out<-stats::integrate(integrand, lower=qnorm(beta), upper=Inf)
    cp<-out$value/(1-beta)
    return(list("n0"=n0, "n1"=n1, "n"=n0+n1, "CP"=round(cp,3)))
  }
  else{
    return(list("error: condition not met. negative n0 returned"))
  }

}
