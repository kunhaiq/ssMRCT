#' @name ssCPsurvPH
#'
#' @title Overall sample size and consistency probability for one MRCT via Method I of MHLW (2007)(conditional version) under random effects model.
#'
#' @description Calculate the overall sample size and consistency probability for one MRCT via Method I of MHLW (2007)(conditional version) under random effects model for survival endpoint with PH assumption.
#'
#' @param alpha one-sided Type I error.
#' @param beta Type II error.
#' @param lambda0r vector of rate parameters in the control group. a \code{R x 1} vector with each component representing the rate parameter in region \eqn{r}.
#' @param hr vector of hazard ratios of treatment group vs control group. a \code{R x 1} vector with each component representing the hr in region \eqn{r}.
#' @param pi threshold ratio in Method I of MHLW (2007)(conditional version). Defaults to 0.5.
#' @param randRatio randomization ratio between the treatment group and control group. Defaults to 1.
#' @param target index for the region of interest. an integer of value referring to the region with the corresponding position in \code{f}.
#' @param f vector of regional fractions. a \code{R x 1} vector with each component representing the regional fraction in region \eqn{r} and the sum must equal 1.
#' @param L scalar of follow-up time in a study assuming fixed study duration design.
#'
#' @details
#' Assume the survival endpoint in region \eqn{r} of group \eqn{k} follows exponential distribution
#' with rate parameter \eqn{\lambda_k^{(r)}}:
#' \eqn{T_k^{(r)} \sim \textrm{exp}(\lambda_k^{(r)}), k=0,1, r=1\ldots R}.
#' Assume the censoring time follows uniform distribution with lower limit \code{unimin} and upper limit \code{unimax}.
#' Assume a study with fixed study duration where each subject is followed for the
#' same length in time (i.e.,\eqn{L} time units) and administrative censoring at the end of study follow-up
#' is the only type of censoring.
#'
#' The overall sample size is calculated based on the following equation,
#' \deqn{
#' \sum_r\left(\tau^2+\frac{\Omega_r}{n_0f_r}\right)^{-1} =\frac{(z_{1-\alpha}+z_{1-\beta})^2}{\delta^2} \quad n_1=\ell n_0.
#' }
#' Since there is no closed form of above equation, \code{ssCPsurvPH} utilizes the \code{\link[stats]{uniroot}} function for numerical solution for \eqn{n_0}.
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
#' Since there is no closed form of above equation, \code{ssCPsurvPH} utilizes the \code{\link[stats]{integrate}} function for numerical integration.
#'
#' Note for survival endpoint under PH assumption,
#' \eqn{\Omega_r=(\ell+1)^2/\{\ell(P_0^{(r)}+\ell P_1^{(r)})\}, P_k^{(r)}=1-e^{-\lambda_k^{(r)L}}}.
#' \eqn{\delta} and \eqn{\tau} are derived based on the naive approach as the mean and standard deviation of the regional hazard ratios.
#'
#'
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
#'            pi=0.5, randRatio=1, target=1, f=rep(1/3,3), L=36)
#'
#'
#' @export
#'
ssCPsurvPH<-function(alpha, beta, lambda0r, hr, pi=0.5, randRatio=1, target, f, L){


  Dr<- -log(hr); delta<-round(mean(Dr),2); tau<-round(sd(Dr),2)

  lambda1r<-lambda0r*hr
  P0<-1-exp(-lambda0r*L)
  P1<-1-exp(-lambda1r*L)
  #Omegars<-4/(P0+P1)
  Omegars<-(randRatio+1)^2/(randRatio*(P0+randRatio*P1))

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

