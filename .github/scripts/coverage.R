# Run from the package root. Both test suites contribute to line coverage.
# Generate reports without RStudio editor integration (also works on CI).
options(covr.rstudio_source_markers = FALSE)
dir.create("coverage", showWarnings = FALSE)
coverage <- covr::package_coverage(type = "tests", quiet = FALSE)
print(coverage)
covr::to_cobertura(coverage, filename = "coverage/cobertura.xml")
saveRDS(coverage, "coverage/coverage.rds")
utils::write.csv(covr::zero_coverage(coverage),
                 "coverage/uncovered-lines.csv", row.names = FALSE)

summary <- c(
  "## Test coverage",
  "",
  sprintf("Overall line coverage: **%.2f%%**.", covr::percent_coverage(coverage)),
  "",
  "Includes `tests/regression.R` and `tests/testthat/`, with no coverage exclusions.",
  "",
  "```text",
  capture.output(print(coverage), type = "message"),
  "```",
  "",
  "Download the test-coverage artifact for uncovered lines and the Cobertura report."
)
writeLines(summary, "coverage/summary.md")
step_summary <- Sys.getenv("GITHUB_STEP_SUMMARY")
if (nzchar(step_summary)) {
  cat(paste(summary, collapse = "\n"), "\n", file = step_summary, append = TRUE)
}
