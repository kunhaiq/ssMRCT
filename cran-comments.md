## Submission

This is the first CRAN submission of ssMRCT (version 0.2.0).
The package calculates regional consistency probabilities and sample sizes
for multi-regional clinical trials under fixed-effects and random-effects
models, including two independent trials and regional treatment-effect estimation.

The license was changed from GPL-3 in the development repository to
MIT + file LICENSE before this first submission. The maintainer email is
kunhaiq@stu.ecnu.edu.cn.

## Local check results

Checked on 2026-10-25 with R 4.5.3, macOS 27.0.1,
aarch64-apple-darwin20, using R CMD check --as-cran on the source tarball.

0 errors | 0 warnings | 1 note

* Unable to verify the current time for the future-file-timestamp check.

Installation, examples, both test suites (including 40 testthat assertions),
README/NEWS checks, PDF manual compilation, HTML validation and HTML math
rendering passed. Pandoc 3.10, HTML Tidy 5.8.0 and V8 were available.

The initial online check completed and reported the expected new-submission
NOTE plus DOI URLs that should use Rd DOI markup. All such help-page URLs
were converted to DOI markup and the manual was regenerated from roxygen2.
The online recheck then stalled during incoming feasibility and was stopped.
The completed final local check used _R_CHECK_CRAN_INCOMING_REMOTE_=false;
remote incoming checks for the corrected tarball still need confirmation.
No code, test or documentation checks were disabled.

## GitHub Actions results

All five R CMD check jobs passed for commit
37cfbad6396baecce5e8175a6b2c3d0a506b93b9 on 2026-10-20:

* Windows, R-release
* macOS, R-release
* Linux, R-release (including PDF manual compilation)
* Linux, R-devel
* Linux, R-oldrel-1

Run: https://github.com/kunhaiq/ssMRCT/actions/runs/39394088679

The coverage workflow also passed:
https://github.com/kunhaiq/ssMRCT/actions/runs/39394088573

These results precede the local DOI-markup correction made on 2026-10-25.
The change replaces ordinary DOI URLs with Rd DOI markup; calculation code and
tests are unchanged. Push the final documentation changes and verify the new
R CMD check run before uploading the source tarball.
