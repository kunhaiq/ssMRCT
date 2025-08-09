#' @name ssCPsurvNPH
#'
#' @title Overall sample size and consistency probability for one MRCT via Japan's criterion I (conditional version) under random effects model.
#'
#' @description Calculate the overall sample size and consistency probability for one MRCT via Japan's criterion I (conditional version) under random effects model for the survival endpoint with nPH assumption.
#'
#' @param alpha Type I error for one-sided test.
#' @param beta Type II error.
#' @param lambda0a vector rate parameter of the piece-wise exponential distribution on or before the breakpoint in the control group. a \code{R x 1} vector with each component representing the rate parameter in region \eqn{r}.
#' @param lambda0b vector rate parameter of the piece-wise exponential distribution after the breakpoint in the control group. a \code{R x 1} vector with each component representing the rate parameter in region \eqn{r}.
#' @param lambda1a vector rate parameter of the piece-wise exponential distribution on or before the breakpoint in the treatment group. a \code{R x 1} vector with each component representing the rate parameter in region \eqn{r}.
#' @param lambda1b vector rate parameter of the piece-wise exponential distribution before the breakpoint in the treatment group. a \code{R x 1} vector with each component representing the rate parameter in region \eqn{r}.
#' @param breakpoint scalar of the time at which the hazard changes. Assumed universal breakpoint value across groups and regions.
#' @param shape0 vector of the shape parameter of the weibull distribution in the control group (shape=1 reduces to exponential). a \code{R x 1} vector with each component representing the shape parameter in region \eqn{r}.
#' @param scale0 vector of the scale parameter of the weibull distribution in the control group (scale=1/rate). a \code{R x 1} vector with each component representing the scale parameter in region \eqn{r}.
#' @param shape1 vector of the shape parameter of the weibull distribution in the treatment group(shape=1 reduces to exponential). a \code{R x 1} vector with each component representing the shape parameter in region \eqn{r}.
#' @param scale1 vector of the scale  parameter of the weibull distribution in the treatment group (scale=1/rate). a \code{R x 1} vector with each component representing the scale parameter in region \eqn{r}.
#' @param PW logical; if TURE (default), survival time follow piece-wise exponential distribution. Otherwise Weibull distribution
#' @param unimin vector of the lower limits of the uniform distribution. a \code{R x 1} vector with each component representing the lower limit in region \eqn{r}.
#' @param unimax vector of the upper limits of the uniform distribution. a \code{R x 1} vector with each component representing the upper parameter in region \eqn{r}.
#' @param pi threshold ratio in Japan's criterion I (conditional version). Defaults to 0.5.
#' @param randRatio randomization ratio between the treatment group and control group. Defaults to 1.
#' @param target index for the region of interest. an integer of value referring to the region with the corresponding position in \code{f}.
#' @param f vector of regional fractions. a \code{R x 1} vector with each component representing the regional fraction in region \eqn{r} and the sum must equal 1.
#' @param eta a scalar of truncated time in RMST. Assumed universal truncated time across groups and regions.
#'
#' @details
#' Assume the survival endpoint follows either piece-wise exponential distribution with hazard function (i.e., rate function)
#' \eqn{\lambda_k^{(r)}I(0<t\le \psi )+\gamma_k^{(r)}I(t> \psi)}, where \eqn{\psi}(\code{breakpoint}) is the change point
#' and \eqn{\lambda_k^{(r)}} (\code{lambda0a}, \code{lambda0b}), \eqn{\gamma_k^{(r)}} (\code{lambda1a}, \code{lambda1b})
#' represent the hazard on or before and after the change point \eqn{\psi} in group \eqn{k} of region \eqn{r}.
#' or Weibull distribution (i.e., Weibull(\eqn{\nu_k^{(r)}, \theta_k^{(r)}})),
#' with shape parameter \eqn{\nu_k^{(r)}} (\code{shape0},\code{shape1})
#' and scale parameter (i.e., 1/rate parameter) \eqn{\theta_k^{(r)}} (\code{scale0}, \code{scale1}) for group \eqn{k} region \eqn{r}.
#' Assume the censoring time follows uniform distribution with lower limit \code{unimin} and upper limit \code{unimax}.
#'
#'
#' The overall sample size is calculated based on the below equation,
#' \deqn{
#' \sum_r\left(\tau^2+\frac{\Omega_r}{n_0f_r}\right)^{-1} =\frac{(z_{1-\alpha}+z_{1-\beta})^2}{\delta^2} \quad n_1=\ell n_0.
#' }
#' Since there is no closed form of above equation, \code{ssCPsurvPH} utilizes the \code{\link[stats]{uniroot}} function for numerical solution for \eqn{n_0}.
#' Then \eqn{n = n_0+n_1}, where both of \eqn{n_0} and \eqn{n_1} should be integers and hence \eqn{n_0}.
#'
#' The consistency probability via Japan's criterion I (conditional version) under random effects model,
#' \eqn{\mathrm{Pr}\left(\widetilde{D}_{r}\ge \pi \widetilde{D} | T>z_{1-\alpha}\right)},
#' is approximately
#' \deqn{
#'  \frac{1}{1-\beta}\int_{u>z_{\beta}}\Phi\left(\frac{(1-\pi) \left(u+z_{1-\alpha}+z_{1-\beta}\right) }{\sqrt{\rho_r^{-1}-1}}\right) d\Phi(u)
#' }
#' where \eqn{\rho_r^{-1}-1=\frac{h_r}{h_r+1}\sum_{j \ne r}^R\frac{h_j}{h_j+1},
#' w=\frac{1}{\tau^2}\sum_{j=1}^R\frac{h_j}{h_j+1}, h_r=\tau^2/\sigma^{2(r)}, \sigma^{2(r)}=\Omega_r/(n_0f_r)}.
#' Since there is no closed form of above equation, \code{ssCPsurvPH} utilizes the \code{\link[stats]{integrate}} function for numerical integration.
#'
#' For survival endpoint under non-PH assumption, \eqn{\Omega_r=\ell^{-1}\sigma_1^{2(r)}+\sigma_1^{2(r)}}.
#' \eqn{\delta} and \eqn{\tau} are derived based on the naive approach as the mean and standard deviation of the regional RMST between group difference.
#'
#'
#' @returns A list containing the following components:
#' \describe{
#'   \item{\code{Dr}}{regional effect}
#'   \item{\code{sigma2.0r}}{asymptotic vairance of RMST estimate of the control group in region r}
#'   \item{\code{sigma2.1r}}{asymptotic vairance of RMST estimate of the treatment group in region r}
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
#' ssCPsurvNPH(alpha=0.025, beta=0.2,
#'             lambda0a=c(0.07,0.07,0.07,0.07),lambda0b=c(0.03,0.04,0.05,0.06),
#'             lambda1a=c(0.02,0.03,0.04,0.05),lambda1b=c(0.03,0.04,0.05,0.06),
#'             breakpoint=10,
#'             unimin=rep(0,4),unimax=rep(80*3,4),
#'             pi=0.5, randRatio=1, target=1, f=rep(1/4,4), eta=80)
#'
#' @export
#'
ssCPsurvNPH<-function(alpha, beta, lambda0a, lambda0b, lambda1a, lambda1b, breakpoint,
                      shape0,scale0,shape1,scale1, PW=T,
                      unimin, unimax, pi=0.5, randRatio=1, target, f, eta){

  rmstTruePW<-function(lambda0a,lambda0b,lambda1a,lambda1b,unimin, unimax, eta, breakpoint){

    ### different in rmst based on the true curve (general weibull distribution)
    integrand<-function(t){
      rmst0<-exp(-lambda0a*t)*(t<=breakpoint)+exp(-lambda0a*breakpoint-lambda0b*(t-breakpoint))*(t>breakpoint);
      rmst1<-exp(-lambda1a*t)*(t<=breakpoint)+exp(-lambda1a*breakpoint-lambda1b*(t-breakpoint))*(t>breakpoint);
      return(rmst1-rmst0)
      #return(exp(-(t/scale1)^shape1)-exp(-(t/scale0)^shape0))
    }
    outw<-integrate(integrand, lower=0, upper=eta)
    d.true<-outw$value


    ###  sigma2up based on true curve (general weibull)
    integrand0<-function(u){  # survival function of group 0
      exp(-lambda0a*u)*(u<=breakpoint)+exp(-lambda0a*breakpoint-lambda0b*(u-breakpoint))*(u>breakpoint)
    }
    inner0<-function(t){

      if(t<=breakpoint){
        (integrate(integrand0, lower=t,upper=eta)$value)^2/(exp(-lambda0a*t)*(unimax-t)/(unimax-unimin))*lambda0a
      }
      else if (t>breakpoint){
        (integrate(integrand0, lower=t,upper=eta)$value)^2/(exp(-lambda0a*breakpoint-lambda0b*(t-breakpoint))*(unimax-t)/(unimax-unimin))*lambda0b
      }

      # (integrate(integrand0, lower=t,upper=eta)$value)^2/((exp(-lambda0a*t)*(t<=breakpoint)+exp(-lambda0a*breakpoint-lambda0b*(t-breakpoint))*(t>breakpoint))*(unimax-t)/(unimax-unimin))*(lambda0a*(t<=breakpoint)+lambda0b*(t>breakpoint))

    }
    result0<-integrate(function(t)sapply(t,inner0),lower=0,upper=eta)    # group 0 sigma2up of RMST


    integrand1<-function(u){ # survival function of group 1
      exp(-lambda1a*u)*(u<=breakpoint)+exp(-lambda1a*breakpoint-lambda1b*(u-breakpoint))*(u>breakpoint)
    }
    inner1<-function(t){
      if(t<=breakpoint){
        (integrate(integrand1, lower=t,upper=eta)$value)^2/(exp(-lambda1a*t)*(unimax-t)/(unimax-unimin))*lambda1a
      }
      else if (t>breakpoint){
        (integrate(integrand1, lower=t,upper=eta)$value)^2/(exp(-lambda1a*breakpoint-lambda1b*(t-breakpoint))*(unimax-t)/(unimax-unimin))*lambda1b
      }

      # (integrate(integrand1, lower=t,upper=eta)$value)^2/((exp(-lambda1a*t)*(t<=breakpoint)+exp(-lambda1a*breakpoint-lambda1b*(t-breakpoint))*(t>breakpoint))*(unimax-t)/(unimax-unimin))*(lambda1a*(t<=breakpoint)+lambda1b*(t>breakpoint))
    }
    result1<-integrate(function(t)sapply(t,inner1),lower=0,upper=eta)   # group 1 sigma2up of RMST

    #sigmaup.true<-sqrt(result0$value/n0+result1$value/n1)
    return(data.frame("dtrue"=d.true, "sigma2up.1"=result1$value,"sigma2up.0"=result0$value))
  }

  rmstTrue1<-function(shape0,scale0,shape1,scale1,unimin, unimax, eta){

    ### validate different in rmst based on the true curve (general weibull distribution)
    integrand<-function(t){
      return(exp(-(t/scale1)^shape1)-exp(-(t/scale0)^shape0))
    }
    outw<-integrate(integrand, lower=0, upper=eta)
    d.true<-outw$value


    ### validate sigma2up based on true curve (general weibull)
    integrand0<-function(u){
      exp(-(u/scale0)^shape0)
    }
    inner0<-function(t){
      (integrate(integrand0, lower=t,upper=eta)$value)^2/(exp(-(t/scale0)^shape0)*(unimax-t)/(unimax-unimin))*shape0/scale0*(t/scale0)^(shape0-1)
    }
    result0<-integrate(function(t)sapply(t,inner0),lower=0,upper=eta)


    integrand1<-function(u){
      exp(-(u/scale1)^shape1)
    }
    inner1<-function(t){
      (integrate(integrand1, lower=t,upper=eta)$value)^2/(exp(-(t/scale1)^shape1)*(unimax-t)/(unimax-unimin))*shape1/scale1*(t/scale1)^(shape1-1)
    }
    result1<-integrate(function(t)sapply(t,inner1),lower=0,upper=eta)

    #sigmaup.true<-sqrt(result0$value/n0+result1$value/n1)
    return(data.frame("dtrue"=d.true, "sigma2up.1"=result1$value,"sigma2up.0"=result0$value))
  }


  regtrue<-NULL; delta.r<-NULL; Omegars<-NULL;sigma2.1r<-NULL;  sigma2.0r<-NULL;
  if(PW==T){  ### PW
    for (r in 1:length(f)){
      regtrue[[r]]<-rmstTruePW(lambda0a=lambda0a[r],lambda0b=lambda0b[r],lambda1a=lambda1a[r],lambda1b=lambda1b[r],
                                 eta=eta,breakpoint=breakpoint,unimin=unimin[r],unimax=unimax[r])

      delta.r[r]<-regtrue[[r]]$dtrue
      sigma2.1r[r]<-regtrue[[r]]$sigma2up.1;
      sigma2.0r[r]<-regtrue[[r]]$sigma2up.0
      #Omegars[r]<-sum(regtrue[[r]][2:3])
      Omegars[r]<-randRatio^(-1)*sigma2.1r[r]+sigma2.0r[r]
    }

  }### end PW ##########

  else{
    #### weibull
    for (r in 1:length(f)){
      regtrue[[r]]<-rmstTrue1(shape0=shape0[r],scale0=scale0[r],shape1=shape1[r],scale1=scale1[r],
                                unimin=unimin[r],unimax=unimax[r],eta=eta)

      delta.r[r]<-regtrue[[r]]$dtrue
      sigma2.1r[r]<-regtrue[[r]]$sigma2up.1;
      sigma2.0r[r]<-regtrue[[r]]$sigma2up.0
      #Omegars[r]<-sum(regtrue[[r]][2:3])
      Omegars[r]<-randRatio^(-1)*sigma2.1r[r]+sigma2.0r[r]
    }
  }### end weibull #####


  delta<-mean(delta.r);
  tau<-sd(delta.r);

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
    return(list("Dr"=delta.r, "sigma2.1r"=sigma2.1r,"sigma2.0r"=sigma2.0r,"n0"=n0, "n1"=n1, "n"=n0+n1, "CP"=round(cp,3)))
  }
  else{
    return(list("error: condition not met. negative n0 returned"))
  }
}



