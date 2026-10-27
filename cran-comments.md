## Resubmission

This is a resubmission of ssMRCT 0.2.0, which has not yet been published on CRAN.

In response to the review dated 2026-10-26, the invalid relative file URI
`LICENSE.md` in README.md has been replaced with the full HTTPS URL:
https://github.com/kunhaiq/ssMRCT/blob/main/LICENSE.md

The same correction was made in README.Rmd. LICENSE.md remains excluded from
the source package; the package's MIT license declaration and LICENSE file
are unchanged. No calculation code or tests were changed.

## Verification of this correction

On 2026-10-27, the source tarball was rebuilt and extracted. With R 4.6.1 and
_R_CHECK_CRAN_INCOMING_CHECK_FILE_URIS_=true, the local incoming checker
reproduced the reported error when the old README link was restored, and
reported no invalid file URIs with the corrected README. The new GitHub URL
returned HTTP 200. This was a targeted link check, not a new full package check.

The same file-URI check is now explicitly enabled in the GitHub Actions workflow.
The earlier full package checks are recorded below; rerun CI after pushing this
README correction before resubmission.

## Package

The package calculates regional consistency probabilities and sample sizes
for multi-regional clinical trials under fixed-effects and random-effects
models, including two independent trials and regional treatment-effect estimation.

The license was changed from GPL-3 in the development repository to
MIT + file LICENSE before this first submission. The maintainer email is
kunhaiq@stu.ecnu.edu.cn.

## Previous local check results (before the README correction)

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
remote incoming checks were not verified in that completed local run.
The subsequent GitHub checks below cover the corrected package.
No code, test or documentation checks were disabled.

## Previous GitHub Actions results (before the README correction)

All five R CMD check jobs passed for commit
152900a33237c6080da954401abfecbbdcafadcd, verified on 2026-10-26:

* Windows, R-release
* macOS, R-release
* Linux, R-release (including PDF manual compilation)
* Linux, R-devel
* Linux, R-oldrel-1

Run: https://github.com/kunhaiq/ssMRCT/actions/runs/40272096589

The R-devel job used R CMD check with --as-cran --no-manual.
The Linux R-release job used --as-cran including PDF manual compilation.
These checks cover the DOI-markup corrections made on 2026-10-25.

The coverage workflow also passed:
https://github.com/kunhaiq/ssMRCT/actions/runs/40272096585
