<!-- README.md is generated from README.Rmd. Please edit that file -->



# ssMRCT

[![R-CMD-check](https://github.com/kunhaiq/ssMRCT/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kunhaiq/ssMRCT/actions/workflows/R-CMD-check.yaml)

**ssMRCT** calculates consistency probabilities (CP), overall sample sizes,
and regional recruitment fractions for multi-regional clinical trials (MRCTs).
It supports fixed-effects and random-effects designs, two independent trials,
and covariate-adjusted regional treatment-effect estimation.

CP is a probability of **regional consistency conditional on overall trial
significance**. It is different from overall power and confidence-interval coverage.
For two trials, the conditioning event is that **both trials are significant**.

## Installation

``` r
# install.packages("pak")
pak::pak("kunhaiq/ssMRCT")
```

## Choose a function

| Task | Functions |
|:-----|:----------|
| Fixed-effects Method I: preserve a fraction of the overall effect | `conProb()`, `conProb2()` |
| Find a regional fraction for a target Method I CP | `regFrac()`, `regFrac2()` |
| Fixed-effects Method II: all regional effects are nonnegative | `conProbII()`, `conProbII2()` |
| Random-effects sample size and Method I CP | `ssCPnorm()`, `ssCPbinary()`, `ssCPsurvPH()`, `ssCPsurvNPH()` |
| Random-effects CP lower bound | `CPlbd()` |
| Covariate-adjusted target-region ATE | `regATE_dr()` |

The suffix `2` denotes two independent trials. `regFrac2()` requires a symmetric
trial design; it does not optimize two different regional fractions.
Random-effects CP uses a **shrinkage regional estimator**, so it has a different
interpretation from fixed-effects CP based on the raw regional effect estimate.

## One trial: CP and regional fraction

Assume a one-sided significance level of 0.05, power of 0.80, treatment effect
of 1, standard deviation of 4 in both arms, and equal randomization. Method I
requires the regional estimate to preserve at least `pi = 0.5` of the overall estimate.


``` r
library(ssMRCT)
conProb(alpha = 0.05, power = 0.8, pi = 0.5,
        rF = 0.271, d = 1, sigmaTrt = 4)
#> $CP
#> [1] 0.8000581
#>
#> $N
#> [1] 396

# Find the regional fraction needed for CP = 0.80.
design <- regFrac(alpha = 0.05, power = 0.8, pi = 0.5, CP = 0.8,
                  d = 1, sigmaTrt = 4)
design
#> $rF
#> [1] 0.2708974
#>
#> $N
#> [1] 396
```

`N` is the recommended overall sample size; `rF` is a continuous regional
fraction. A recruitment plan must allocate whole participants to each regional
arm. Multiplying `N` by `rF` alone does not establish CP for that integer allocation.

## Two trials: general CP and symmetric fraction selection

The CP calculation allows different trial powers and regional fractions.
Here study 2 has higher power; unspecified effect and variance parameters
are inherited from study 1.


``` r
conProb2(alpha = 0.05, power1 = 0.8, power2 = 0.9,
         rF1 = 0.15, rF2 = 0.20, d1 = 1, sigmaTrt1 = 4)
#> $CP
#> [1] 0.8309535
#>
#> $N1
#> [1] 396
#>
#> $N2
#> [1] 550
```

For two symmetric trials, find a common fraction and verify its CP:


``` r
fraction <- regFrac2(alpha = 0.05, power = 0.8, pi = 0.5, CP = 0.8)
fraction
#> [1] 0.1530893
conProb2(alpha = 0.05, power1 = 0.8, pi = 0.5,
         rF1 = fraction, d1 = 1, sigmaTrt1 = 4)
#> $CP
#> [1] 0.8
#>
#> $N1
#> [1] 396
#>
#> $N2
#> [1] 396
```

Both trials must share the effect, arm-specific variances, randomization ratio,
power and regional fraction for this inverse calculation. For unequal
randomization ratios in `conProb2()` or continuous `conProbII2()`, specify the
control means (`muCtrl1`, `muCtrl2`) because separate arm pooling can make them
relevant to the pooled effect.

## All-region consistency with a binary endpoint

Method II requires every regional treatment-minus-control estimate to be
nonnegative. This example uses three regions, response probabilities of 0.70
and 0.55, and 2:1 treatment-to-control randomization.


``` r
set.seed(42)
binary <- conProbII(alpha = 0.025, power = 0.8, rF = c(0.2, 0.4, 0.4),
                    pTrt = 0.70, pCtrl = 0.55, randRatio = 2,
                    responseType = "binary", B = 50000)
binary[c("CP", "N", "nTrtReg", "nCtrlReg", "achievedPower", "MCSE")]
#> $CP
#> [1] 0.9151825
#>
#> $N
#> [1] 373
#>
#> $nTrtReg
#> [1] 50 99 99
#>
#> $nCtrlReg
#> [1] 25 50 50
#>
#> $achievedPower
#> [1] 0.80856
#>
#> $MCSE
#> [1] 0.001385656
```

The simulation uses the returned integer regional arm counts. CP is the
proportion of significant trials that also meet regional consistency;
zero regional effects count as consistent. `achievedPower` is the simulated
probability of overall significance, and `MCSE` describes simulation uncertainty
in CP. For `conProbII2()`, `achievedPower` instead refers to both trials succeeding.

## Random-effects design

For a normal endpoint, supply the mean regional effect `delta`, between-region
standard deviation `tau`, outcome standard deviations, and regional fractions.


``` r
ssCPnorm(alpha = 0.025, beta = 0.2, delta = 0.25, tau = 0.1,
         sigma1r = 1, sigma0r = 1, randRatio = 1,
         pi = 0.5, target = 1, f = c(0.2, 0.3, 0.5))
#> $n0
#> [1] 459
#>
#> $n1
#> [1] 459
#>
#> $n
#> [1] 918
#>
#> $CP
#> [1] 0.9934161
```

`n0`, `n1` and `n` are **overall** control, treatment and total counts.
`CP` refers to the region indexed by `target`. The function sizes the trial for
specified overall power and then evaluates regional CP; it does not solve for
a user-specified CP target. A finite sample-size solution requires
`tau * (qnorm(1 - alpha) + qnorm(1 - beta)) / delta < sqrt(length(f))`.

`ssCPbinary()` uses regional response probabilities and a normal approximation.
`ssCPsurvPH()` uses exponential survival with fixed follow-up; `ssCPsurvNPH()`
uses RMST differences with piecewise exponential or Weibull survival and uniform
censoring. These functions derive `delta` and `tau` from the unweighted mean and
sample standard deviation of the supplied regional effects.

`CPlbd()` returns a CP lower bound and **regional** counts at its equality
condition, not an overall allocation. Its `fr` argument specifies output length;
it does not divide the counts by regional fractions. Counts are `NA` when the
finite equality condition does not apply.

## Calculation conventions

- `alpha` is one-sided. Fixed-effects functions take `power`; random-effects
  functions take `beta = 1 - power`.
- `randRatio` is treatment/control. `sigma` parameters are standard deviations,
  not variances. Higher effects indicate benefit; PH effects are `-log(hr)`.
- Method I uses `pi` for the fraction of the overall effect to preserve.
  In `regATE_dr()`, `pi` instead denotes target-region membership probability.
- Analytical CP uses the continuous sample-size design. Returned integer counts
  are recommendations; CP is not recalculated after regional rounding.
  Binary Method II simulation uses the actual returned integer allocation.
- Set a seed for binary simulations and continuous Method II numerical
  integration. The latter uses randomized multivariate normal integration.

## Regional ATE estimation

`regATE_dr()` accepts participant-level data with `y`, `trt`, `region`, and
numeric covariates. It fits outcome and probability models to borrow information
for a target-region ATE. With `region_wise = TRUE`, it screens source regions
using arm-specific ANCOVA comparisons before refitting.

Borrowing requires conditional outcome-mean exchangeability across retained
regions. The reported influence-function interval does not include uncertainty
from source-region selection. See `?regATE_dr` for the data format and model options.

## Documentation and references

Use `help(package = "ssMRCT")` for the function index and `?conProb` (or another
function name) for arguments, formulas, return values and examples.

- Qing et al. (2025), regional consistency under two MRCTs:
  [doi:10.1002/sim.70306](https://doi.org/10.1002/sim.70306).
- Ren and Xu (2026), random-effects consistency and regional sample size:
  [doi:10.1002/pst.70090](https://doi.org/10.1002/pst.70090).

See [NEWS.md](NEWS.md) for changes that affect calculations and returned values.

## License

MIT license. See [LICENSE.md](LICENSE.md) for the full text.

## Automated checks

GitHub Actions runs package checks on every push and pull request, and can also
be started from the Actions tab. The matrix covers Windows and macOS with
R-release, plus Linux with R-release, R-devel and the previous R release.
Errors and warnings fail the check; notes remain in the logs for review.
Check reports are saved as workflow artifacts, including successful runs.

The matrix checks examples, regression tests and help files using `--as-cran`.
The Linux R-release job also compiles the PDF manual. All check logs and notes
must be reviewed before a CRAN submission. The workflow does not publish the package to CRAN.
