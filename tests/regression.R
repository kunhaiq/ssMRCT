library(ssMRCT)
near <- function(x, y, tol=1e-6) stopifnot(isTRUE(all.equal(x,y,tolerance=tol,check.attributes=FALSE)))
fails <- function(expr, pattern) {
  msg <- tryCatch({force(expr); NA_character_}, error=function(e) conditionMessage(e))
  stopifnot(!is.na(msg), grepl(pattern,msg))
}
# Published benchmarks and inverse-design consistency.
x <- conProb(.05,.8,rF=.271,d=1,sigmaTrt=4)
near(x$CP, .8000581); stopifnot(x$N == 396)
f <- regFrac(.05,.8)$rF
near(conProb(.05,.8,rF=f,d=1,sigmaTrt=4)$CP,.8)
f2 <- regFrac2(.05,.8)
for (d in c(1,4,8)) near(conProb2(.05,.8,rF1=f2,d1=d,sigmaTrt1=4)$CP,.8)
near(conProb(.05,.8,rF=1,d=1,sigmaTrt=4)$CP,1)
near(conProb2(.05,.8,rF1=1,d1=1,sigmaTrt1=4)$CP,1)
stopifnot(regFrac2(.05,.8,CP=1)==1)
fails(conProb(.05,.8,rF=.2,d=-1,sigmaTrt=4),'d')
fails(regFrac(.05,.8,CP=.4),'CP')
set.seed(21)
near(conProbII(.05,.8,rF=rep(1/3,3),d=1,sigmaTrt=4)$CP,.89072,1e-4)
near(conProbII2(.05,.8,rF1=rep(1/3,3),d1=1,sigmaTrt1=4)$CP,.9836,2e-4)
# Baselines matter when arm pooling weights differ; common translation does not.
args <- list(alpha=.05,power1=.8,rF1=c(.2,.8),rF2=c(.4,.6),d1=1,
             sigmaTrt1=4,randRatio1=1,randRatio2=2)
set.seed(3); baseline <- do.call(conProbII2,args)$CP
set.seed(3); shifted <- do.call(conProbII2,c(args,list(muCtrl1=10,muCtrl2=10)))$CP
near(baseline,shifted,1e-5)
set.seed(3); changed <- do.call(conProbII2,c(args,list(muCtrl1=0,muCtrl2=-20)))$CP
stopifnot(abs(baseline-changed)>.01)
# With one region, overall significance implies regional consistency exactly.
for (p in list(c(.7,.5),c(.9,.1))) {
 set.seed(789)
 x <- conProbII(.05,.8,rF=1,pTrt=p[1],pCtrl=p[2],responseType='binary',B=2000)
 stopifnot(x$CP==1,x$MCSE==0,x$successes==x$consistentSuccesses)
 x <- conProbII2(.05,.8,rF1=1,pTrt1=p[1],pCtrl1=p[2],responseType='binary',B=2000)
 stopifnot(x$CP==1,is.finite(x$achievedPower))
}
# Independent exact enumeration at a small two-region design; ties have positive mass.
f <- c(.5,.5);pT <- .8;pC <- .2
set.seed(91)
x <- conProbII(.05,.8,rF=f,pTrt=pT,pCtrl=pC,responseType='binary',B=100000)
nT <- x$nTrtReg;nC <- x$nCtrlReg
states <- expand.grid(t1=0:nT[1],t2=0:nT[2],c1=0:nC[1],c2=0:nC[2])
pt <- (states$t1+states$t2)/sum(nT);pc <- (states$c1+states$c2)/sum(nC)
se <- sqrt(pt*(1-pt)/sum(nT)+pc*(1-pc)/sum(nC))
success <- pt-pc > qnorm(.95)*se
consistent <- states$t1*nC[1]>=states$c1*nT[1] & states$t2*nC[2]>=states$c2*nT[2]
prob <- dbinom(states$t1,nT[1],pT)*dbinom(states$t2,nT[2],pT)*
        dbinom(states$c1,nC[1],pC)*dbinom(states$c2,nC[2],pC)
exact <- sum(prob[success & consistent])/sum(prob[success])
stopifnot(abs(x$CP-exact) < 5*x$MCSE+.001)
set.seed(432)
x <- conProbII(.05,.8,rF=c(.2,.4,.4),pTrt=.8,pCtrl=.5,randRatio=2,
               responseType='binary',B=50000)
stopifnot(x$nTrt==48,x$nCtrl==25,x$N==73,abs(x$CP-.8912)<.012)
fails(conProbII(.05,.8,rF=c(.2,.5),pTrt=.8,pCtrl=.5,responseType='binary'),'sum')
fails(conProbII(.05,.8,rF=1,pTrt=.8,pCtrl=.5,responseType='binary',B=1.5),'integer')
# Closed-form homogeneous random-effects sample size, including roots outside the old interval.
a <- (qnorm(.975)+qnorm(.8))^2
for (par in list(c(.25,.1),c(.01,.001),c(10,.1))) {
 x <- ssCPnorm(.025,.2,par[1],par[2],1,1,target=1,f=rep(1/3,3))
 expected <- max(1,ceiling(6/(3*par[1]^2/a-par[2]^2)))
 stopifnot(x$n0==expected,x$n1==expected,is.finite(x$CP))
}
fails(ssCPnorm(.025,.2,.1,1,1,1,target=1,f=rep(1/3,3)),'No finite')
x <- CPlbd(.05,.2,.4,.25,1,1,c(.1,.5))
stopifnot(all(x$n0r==196),all(x$nr==392))
x <- CPlbd(.05,.2,1,.25,1,1,c(.1,.5))
stopifnot(is.numeric(x$n0r),all(is.na(x$n0r)))
# Bernoulli moments must reach the common normal design without intermediate rounding.
p0 <- c(.3,.4,.5);p1 <- c(.407,.51,.613)
x <- ssCPbinary(.025,.2,p0,p1,target=2,f=rep(1/3,3))
y <- ssCPnorm(.025,.2,mean(p1-p0),sd(p1-p0),sqrt(p1*(1-p1)),sqrt(p0*(1-p0)),
              target=2,f=rep(1/3,3))
near(unlist(x),unlist(y))
# No censoring before eta: Var(min(T,eta)) from its first two moments.
rmst_moments <- function(lambda, eta) {
 m <- (1-exp(-lambda*eta))/lambda
 second <- 2*(1-exp(-lambda*eta)*(1+lambda*eta))/lambda^2
 second-m^2
}
x <- ssCPsurvNPH(.025,.2,lambda0a=.1,lambda0b=.1,lambda1a=.07,lambda1b=.07,
  breakpoint=5,unimin=20,unimax=30,target=1,f=rep(1/3,3),eta=10)
near(x$sigma2.0r,rep(rmst_moments(.1,10),3));near(x$sigma2.1r,rep(rmst_moments(.07,10),3))
stopifnot(x$n==514)
y <- ssCPsurvNPH(.025,.2,shape0=1,scale0=10,shape1=1,scale1=1/.07,PW=FALSE,
  unimin=20,unimax=30,target=1,f=rep(1/3,3),eta=10)
near(unlist(x),unlist(y))
# ATE scalar/vector agreement and target-only fallback after source exclusion.
set.seed(123)
dat <- data.frame(region=rep(1:3,each=100),trt=rep(c(0,1),150),x=rnorm(300))
dat$y <- dat$x+dat$trt+rnorm(300)
x <- regATE_dr(dat,p1=.5,p0=.5,pi=1/3)
y <- regATE_dr(dat,p1=rep(.5,300),p0=rep(.5,300),pi=rep(1/3,300))
near(x$estimate,y$estimate);near(x$se_eif,y$se_eif)
dat$y <- dat$y+100*(dat$region!=1)
y <- regATE_dr(dat,p1=rep(.5,300),p0=rep(.5,300),pi=rep(1/3,300),region_wise=TRUE)
stopifnot(length(y$region_wise_remaining_regions)==0,is.finite(y$region_wise_estimate))
fails(regATE_dr(dat,p1=c(.5,.5)),'nrow')
# Ren et al. (2026), Corollary 2, equation (15), independent of the sample-size solver.
for (R in c(3,4)) {
 z <- qnorm(.975)+qnorm(.8);delta <- .25;tau <- .1
 ref <- integrate(function(u) pnorm(delta^2*R*.5*(u+z)/(tau^2*sqrt(R-1)*z^2))*dnorm(u),
                  qnorm(.2),Inf)$value/.8
 ans <- ssCPnorm(.025,.2,delta,tau,1,1,target=1,f=rep(1/R,R))
 near(ans$CP,ref,1e-7)
}
# Equation (7): at tau=0 the regional shrinkage estimate equals the overall one.
near(ssCPnorm(.025,.2,.25,0,1,1,target=1,f=c(.1,.3,.6))$CP,1)
# Independent Gaussian simulation of equations (2) and (7), unequal regional variances.
f <- c(.1,.3,.6);Omega <- c(2,3,5);delta <- .25;tau <- .1;z <- qnorm(.975)+qnorm(.8)
n <- uniroot(function(n) sum(1/(tau^2+Omega/(n*f)))-z^2/delta^2,c(1,10000))$root
variance <- Omega/(n*f);w <- 1/(tau^2+variance)
set.seed(802)
b <- 150000
D <- sapply(tau^2+variance,function(v) rnorm(b,delta,sqrt(v)))
overall <- as.vector(D %*% w/sum(w))
regional <- tau^2/(tau^2+variance[1])*D[,1]+variance[1]/(tau^2+variance[1])*overall
success <- overall*sqrt(sum(w))>qnorm(.975)
empirical <- mean((regional >= .5*overall)[success])
ans <- ssCPnorm(.025,.2,delta,tau,sqrt(Omega/2),sqrt(Omega/2),target=1,f=f)
stopifnot(abs(empirical-ans$CP)<.005)
stopifnot(ans$CP >= CPlbd(.025,.2,tau/delta,delta,1,1,.1)$CPlbd-1e-6)
# Partial source removal: row-aligned treatment scores and re-estimated selection scores.
set.seed(903)
dat <- data.frame(region=rep(1:3,each=200),trt=rep(c(0,1),300),x=rnorm(600))
dat$y <- dat$x+dat$trt+rnorm(600)+100*(dat$region==2)
x <- suppressWarnings(regATE_dr(dat,p1=rep(.5,600),p0=rep(.5,600),
                                pi=rep(1/3,600),region_wise=TRUE))
y <- regATE_dr(dat[dat$region!=2,],p1=.5,p0=.5)
stopifnot(identical(as.character(x$region_wise_remaining_regions),'3'))
near(x$region_wise_estimate,y$estimate);near(x$region_wise_se,y$se_eif)
# Full-region Method I is exactly consistent, even near complete preservation.
near(conProb2(.05, .8, pi = .99999999, rF1 = 1, d1 = 1,
              sigmaTrt1 = 4)$CP, 1)
# Very small heterogeneity must retain finite lower-bound equality counts.
x <- CPlbd(.05, .2, 1e-200, .25, 1, 1, .2)
near(x$CPlbd, 1)
stopifnot(x$n0r == ceiling(2 * ((qnorm(.95) + qnorm(.8)) / .25)^2 / 2))
# Large post-breakpoint hazards must not create Inf * 0 in inactive branches.
x <- ssCPsurvNPH(.025, .2, lambda0a = .001, lambda0b = 1000,
  lambda1a = .0007, lambda1b = 1000, breakpoint = 5,
  unimin = 20, unimax = 30, target = 1, f = rep(1/3, 3), eta = 10)
stopifnot(all(is.finite(unlist(x))), all(x$sigma2.0r > 0), x$n > 0)
# Reject ambiguous labels, unsupported settings and unusable covariates early.
fails(regATE_dr(dat, target_region = c(1, 2)), 'single nonmissing')
fails(regATE_dr(dat, config = list(unknown = 1)), 'check_p_threshold')
fails(regATE_dr(dat, compute_se = NA), 'TRUE or FALSE')
bad <- dat; bad$x[1] <- Inf
fails(regATE_dr(bad), 'finite')
bad <- dat; bad$D <- 0
fails(regATE_dr(bad), 'reserved')

# Unequal two-trial design: compare CP with independently generated arm means.
# Generate independent regional arm means, then construct the observed endpoints.
set.seed(20261013)
B <- 200000
f1 <- c(.2,.3,.5); f2 <- c(.3,.4,.3)
simulate_trial <- function(power, effect, st, sc, ratio, baseline, fractions) {
 z <- qnorm(.975) + qnorm(power)
 nc <- (st^2 / ratio + sc^2) * z^2 / effect^2
 nt <- ratio * nc
 yt <- sapply(fractions, function(f) rnorm(B, baseline + effect, st/sqrt(nt*f)))
 yc <- sapply(fractions, function(f) rnorm(B, baseline, sc/sqrt(nc*f)))
 overall <- as.vector((yt-yc) %*% fractions)
 list(yt=yt,yc=yc,nt=nt*fractions,nc=nc*fractions,
      success=overall>qnorm(.975)*sqrt(st^2/nt+sc^2/nc))
}
a <- simulate_trial(.8,1,4,3.5,1,0,f1)
b <- simulate_trial(.9,.8,3.5,3,2,.5,f2)
pooledT <- (sweep(a$yt,2,a$nt,`*`)+sweep(b$yt,2,b$nt,`*`))
pooledC <- (sweep(a$yc,2,a$nc,`*`)+sweep(b$yc,2,b$nc,`*`))
regional <- sweep(pooledT,2,a$nt+b$nt,`/`)-sweep(pooledC,2,a$nc+b$nc,`/`)
overall <- rowSums(pooledT)/sum(a$nt+b$nt)-rowSums(pooledC)/sum(a$nc+b$nc)
success <- a$success & b$success
emp1 <- mean((regional[,1]>=.5*overall)[success])
emp2 <- mean((rowSums(regional>=0)==3)[success])
args <- list(alpha=.025,power1=.8,power2=.9,d1=1,d2=.8,
 sigmaTrt1=4,sigmaCtrl1=3.5,sigmaTrt2=3.5,sigmaCtrl2=3,
 randRatio1=1,randRatio2=2,muCtrl1=0,muCtrl2=.5)
cp1 <- do.call(conProb2,c(args,list(rF1=f1[1],rF2=f2[1])))$CP
cp2 <- do.call(conProbII2,c(args,list(rF1=f1,rF2=f2)))$CP
mcse <- sqrt(c(emp1*(1-emp1),emp2*(1-emp2))/sum(success))
stopifnot(all(abs(c(cp1-emp1,cp2-emp2)) < 5*mcse + 1e-5))
