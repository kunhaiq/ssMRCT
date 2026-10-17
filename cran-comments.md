## Submission

This is the first CRAN submission of ssMRCT (version 0.2.0).
The package calculates regional consistency probabilities and sample sizes
for multi-regional clinical trials under fixed-effects and random-effects
models, including two independent trials.

The license was changed from GPL-3 in the development repository to
MIT + file LICENSE before this first submission. The maintainer email is
kunhaiq@stu.ecnu.edu.cn.

## Local check results

Checked on 2026-10-17 with R 4.5.3, macOS 27.0.1, aarch64-apple-darwin20,
using R CMD check --as-cran on the source tarball.

0 errors | 0 warnings | 1 note

* Unable to verify the current time for the future-file-timestamp check.

Pandoc, HTML Tidy 5.8.0 and V8 were available. Installation, examples,
regression tests, README/NEWS checks, PDF manual compilation, HTML validation
and HTML math rendering passed.

The full online run stalled during CRAN incoming feasibility and was stopped.
The completed local run used _R_CHECK_CRAN_INCOMING_REMOTE_=false; remote
incoming checks have therefore not been fully verified in this run. No code,
test or documentation checks were disabled.

## Automated checks prepared; remote results pending

.github/workflows/R-CMD-check.yaml runs on every push and pull request,
and can be started manually. Its matrix covers:

* Windows, R-release
* macOS, R-release
* Linux, R-release (including PDF manual compilation)
* Linux, R-devel
* Linux, R-oldrel-1

These checks have been configured but have not yet run on GitHub.
Before submitting, push the prepared changes, inspect all five job results
and their notes, and replace this section with the actual run URL and results.
The local R 4.5.3 check does not replace a current R-release/R-devel check.
This remains a preparation draft until those results are recorded.
