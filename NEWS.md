# ssMRCT 0.2.0

* Fix the README license link to use its full GitHub URL, since LICENSE.md
  is not included in the CRAN source package. Enable file-URI checks in CI.

* Use DOI markup in help-page references as required by CRAN incoming checks.

* Disable RStudio editor integration when generating coverage reports so the
  GitHub Actions coverage job does not require rstudioapi.

* Add GitHub Actions line-coverage reports for both test suites, with downloadable
  uncovered-line reports and optional Codecov uploads.

* Add grouped testthat tests for common trial designs, survival moments and
  regional ATE estimation, alongside the existing regression suite.

* Label consistency and design help-page titles by Method I or Method II,
  and remove print-method arguments from the regional ATE help page.

* Use complete Chicago author-date references in the README and help pages,
  including all authors, journal volume and issue, article identifiers and DOIs.

* Expand the regional ATE README guide with participant-level data, model options,
  and executable examples comparing estimation before and after source screening.

* Add GitHub Actions checks on push and pull requests for Windows, macOS and
  Linux, including R-devel and the previous R release. The Linux R-release
  job also compiles the PDF manual; successful check logs are retained.

* Prepare the first CRAN submission: adopt the MIT license, add structured
  author metadata and method references, and update the maintainer email to
  kunhaiq@stu.ecnu.edu.cn.

* Correct binary Method II arm allocation for unequal randomization, use actual
  arm sizes in Wald standard errors, include zero regional effects, and estimate
  conditional CP using observed overall successes. Return simulation diagnostics
  and regional arm sizes. Zero-variance samples are handled without division by zero.
* Accept the default continuous response type in both Method II functions.
* Use continuous design sample sizes consistently inside theoretical CP formulas;
  return rounded overall sample size recommendations separately. `conProb2()` and
  random-effect CP values may change slightly. `regFrac2()` remains a scalar-returning
  solver for the symmetric design, not a general unconstrained optimizer.
* Allow different continuous control means in `conProbII2()` via trailing
  `muCtrl1` and `muCtrl2` arguments, preserving earlier positional arguments.
* Share design validation and random-effect root finding. Search beyond the old
  [1, 100000] bracket, handle small roots, and distinguish infeasibility from invalid
  inputs. Errors now raise conditions rather than returning unnamed message lists.
* Preserve internal precision for Bernoulli and PH moments and CP outputs.
* Correct uniform censoring survival before its lower limit for NPH RMST variance.
  Require the RMST horizon to be strictly below the censoring upper limit.
* Correct `CPlbd()` regional sample sizes (remove the erroneous inverse-fraction
  factor); return numeric NA for unavailable counts. The `fr` argument remains
  accepted for compatibility and regional vector sizing.
* Fix supplied propensity vectors in `regATE_dr()`, subset treatment propensity
  vectors during source selection, re-estimate region propensity when the retained
  population changes, and handle target-only estimation explicitly.
* Add regression tests against published benchmarks, exact binomial probabilities,
  no-censoring RMST moments, and Ren et al. (2026) random-effects identities and simulation.
* Simplify public roxygen documentation, explain parameter conventions and returned
  counts, distinguish continuous-design CP from rounded recommendations, and clarify
  source-screening assumptions and uncertainty in regional ATE estimation.
* Return exact full-region Method I CP in the degenerate two-trial case, and avoid
  numerical overflow in piecewise-hazard RMST integration and small-heterogeneity
  lower-bound equality counts.
* Validate ATE settings, target labels and covariates, and reject undefined nuisance
  predictions and weights with explicit errors.
* Remove the unused cubature dependency.
* Rewrite the README with a function guide and executable design examples; clarify
  CP conditioning, symmetric fraction selection, and overall versus regional counts.
* Verify two-trial Method I and II formulas against independent normal simulations
  with unequal powers, variances, regional fractions and randomization ratios.
