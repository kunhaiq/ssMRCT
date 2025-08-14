
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ssMRCT

<!-- badges: start -->

<!-- badges: end -->

The goal of ssMRCT is to provide a comprehensive suite of functions for
designing and analyzing Multi-Regional Clinical Trials (MRCTs) under
fixed effects model and random effects model for regulatory decision
support and trial optimization in global drug development programs.

## Installation

You can install the development version of ssMRCT from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("kunhaiq/ssMRCT")
```

## Example

This is a basic example which shows you how to solve a common problem.
Suppose the significance level is 0.05 and the desired power is 0.8.
Additionally, the difference in mean of the clinical trial response is 1
and its corresponding variance is 16. We need to compute the minimal
overall sample size and regional sample size, which equals the regional
proportion times the overall sample size.

``` r
library(ssMRCT)
## conProb() gives the consistency probability under fixed effects model
conProb(alpha = 0.05, power = 0.8, rF = 0.271, d = 1, sigmaTrt = 4)
#> $CP
#> [1] 0.8000581
#> 
#> $N
#> [1] 396
## regFrac() gives the optimal regional fraction under fixed effects model
regFrac(alpha = 0.05, power = 0.8, d = 1, sigmaTrt = 4)
#> $rF
#> [1] 0.2708725
#> 
#> $N
#> [1] 396
```

Also we could consider to conduct two MRCTs

``` r
## conProb2() gives the consistency probability under fixed effects model
conProb2(alpha = 0.05, power1 = 0.8, power2 = 0.9, rF1 = 0.141, d1 = 1, sigmaTrt1 = 4) 
#> $CP
#> [1] 0.8002532
#> 
#> $N1
#> [1] 396
#> 
#> $N2
#> [1] 550
## regFrac2() gives the optimal regional fractions under fixed effects model
regFrac2(alpha = 0.05, power1 = 0.8, power2 = 0.9, d1 = 1, sigmaTrt1 = 4)
#> $rF1
#> [1] 0.1407622
#> 
#> $rF2
#> [1] 0.1407622
#> 
#> $N1
#> [1] 396
#> 
#> $N2
#> [1] 550
```

For MRCT under random effects model, one could use

``` r
## ssCPnorm() gives the consistency probability under random effects model
ssCPnorm(alpha=0.025, beta=0.2, delta=0.25, tau=0.1, sigma1r=rep(1,3),
         sigma0r=rep(1,3), randRatio=1, pi=0.5, target=1, f=c(1/3,1/3,1/3))
#> $n0
#> [1] 433
#> 
#> $n1
#> [1] 433
#> 
#> $n
#> [1] 866
#> 
#> $CP
#> [1] 0.99
```
