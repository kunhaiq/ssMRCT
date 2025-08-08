#' @name ssCPsurvPH
#'
#' @title Overall sample size and consistency probability for one MRCT via Japan's criterion I (conditional version) under random effects model.
#'
#' @description Calculate the overall sample size and consistency probability for one MRCT via Japan's criterion I (conditional version) under random effects model for survival endpoint with PH assumption.
#'
#' @param alpha one-sided Type I error.
#' @param beta Type II error. 
#' @param lambda0r vector rate parameter in the exponential distribution in the control group for all regions.
#' @param hr vector of hazard ratios of treatment group vs control group for all regions.
#' @param pi threshold ratio in Japan's criterion I (conditional version). Defaults to 0.5.
#' @param k randomization ratio between the treatment group and control group. Defaults to 1.
#' @param target index for the region of interest. Should be an integer.
#' @param f vector of regional fractions.Specified for all regions with length that equals to the total number of region.
#' @param L scalar of follow-up time in a fixed study duration design.
#'
#' @details
#' The overall sample size is calculated based on the below equation,
#' \deqn{
#' \sum_r\left(\tau^2+\frac{\Omega_r}{n_0f_r}\right)^{-1} =\frac{(z_{1-\alpha}+z_{1-\beta})^2}{\delta^2} \quad n_1=kn_0.
#' }
#' Since there is no closed form of above equation, \code{ssCPsurvPH} utilizes the \code{\link[stats]{uniroot}} function for numerical solution for \eqn{n_0}.
#' Then \eqn{n = n_0+n_1}. Additionally, both of
#' \eqn{n_0} and \eqn{n_1} should be integers and hence \eqn{n_0}.
#'
#' The consistency probability via Japan's criterion I (conditional version) under random effects model,
#' \eqn{\mathrm{Pr}\left(\widetilde{D}_{r}\ge \pi \widetilde{D} | T>z_{1-\alpha}\right)},
#' is approximately
#' \deqn{
#'  \frac{1}{1-\beta}\int_{u>z_{\beta}}\Phi\left(\frac{(1-\pi) \left(u+z_{1-\alpha}+z_{1-\beta}\right) }{\sqrt{\rho_r^{-1}-1}}\right) d\Phi(u)
#' }provided \eqn{\frac{\tau}{\delta}<\frac{\sqrt{R}}{z_{1-\alpha}+z_{1-\beta}}} where \eqn{\rho_r^{-1}=\rho_r^{-1}=1+\frac{h_r}{h_r+1}\sum_{j \ne r}^R\frac{h_j}{h_j+1},
#' \quad w=\frac{1}{\tau^2}\sum_{j=1}^R\frac{h_j}{h_j+1}, \quad h_r=\tau^2/\sigma^{2(r)},
#' \quad \sigma^{2(r)}=\Omega_r/(n_0f_r)}
#' Since there is no closed form of above equation, \code{ssCPsurvPH} utilizes the \code{\link[stats]{integrate}} function for numerical integration.
#'
#' Assume the survival endpoint follows exponential distribution \eqn{T_k^{(r)} \sim \textrm{exp}(\lambda_k^{(r)}), k=0,1; r=1\ldots R}.
#' 1) \eqn{\omega_r=(\ell+1)^2/\{\ell(P_0^{(r)+\ell P_1^{(r)}})\}, \quad R_k^{(r)=1-e^{-\lambda_k^{(r)L}}}}
#' 2) \eqn{\delta} and \eqn{\tau} are derived based on the naive approach as the mean and standard deviation of the regional hazard ratios.
#'
#' Assume the censoring time follows uniform distribution with lower limit \code{unimin} and upper limit \code{unimax}
#'
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
#' ssCPsurvPH(alpha=0.025, beta=0.2, lambda0r=rep(0.05,3), hr=c(0.7,0.6,0.4), 
#'             pi=0.5, k=1, target=1, f=rep(1/3,3), L=36)
#'
#' @export
#'
ssCPsurvPH<-function(alpha, beta, lambda0r, hr, pi, k, target, f, L){


  Dr<- -log(hr); delta<-round(mean(Dr),2); tau<-round(sd(Dr),2)

  lambda1r<-lambda0r*hr
  P0<-1-exp(-lambda0r*L)
  P1<-1-exp(-lambda1r*L)
  #Omegars<-4/(P0+P1)
  Omegars<-(k+1)^2/(k*(P0+k*P1))

  R<-length(f);
  tau_delta_r<-tau/delta
  condition<-sqrt(R)/(qnorm(1-alpha)+qnorm(1-beta))  ## upper value of tau_delta_r

  if(abs(sum(f)-1) > 10^(-10)){return(list("error: regional fractions do not sum up to 1"))}
  else if (tau_delta_r<condition){ #### add a check to ensure positive n0 to fit into the context

    #### calculate overall sample size n0
    fp<-function(n0){
      rhs <- (qnorm(1-alpha)+qnorm(1-beta))^2/delta^2
      sum((tau^2+Omegars/(n0*f))^(-1))-rhs
    }
    n0<-ceiling(stats::uniroot(fp,c(1,100000), tol = 0.0001)$root) 	# solve eq.
    n1<-ceiling(n0*k)
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
    return(list("n0"=n0, "n1"=n1, "n"=n0+n1, "CP"=round(cp,5)*100))
  }
  else{
    return(list("error: condition not met. negative n0 returned"))
  }
}

 