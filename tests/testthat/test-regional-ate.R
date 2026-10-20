# Balanced randomized data, with a shared effect and one baseline covariate.
make_ate_data <- function() {
  set.seed(42)
  dat <- data.frame(
    region = rep(1:3, each = 100),
    trt = unlist(lapply(1:3, function(i) sample(rep(c(0, 1), each = 50)))),
    x = rnorm(300)
  )
  dat$y <- dat$x + dat$trt + rnorm(300)
  dat
}

test_that("ATE and uncertainty transform correctly when outcome units change", {
  dat <- make_ate_data()
  a <- regATE_dr(dat, p1 = 0.5, p0 = 0.5)
  dat$y <- 10 + 2 * dat$y
  b <- regATE_dr(dat, p1 = 0.5, p0 = 0.5)
  expect_equal(b$estimate, 2 * a$estimate, tolerance = 1e-8)
  expect_equal(b$se_eif, 2 * a$se_eif, tolerance = 1e-8)
  expect_equal(b$ci, 2 * a$ci, tolerance = 1e-8)
})

test_that("source screening matches a separate fit to the retained regions", {
  dat <- make_ate_data()
  dat$y <- dat$y + 1.5 * (dat$region == 3) * dat$trt
  screened <- regATE_dr(dat, p1 = 0.5, p0 = 0.5, region_wise = TRUE)
  retained <- regATE_dr(dat[dat$region != 3, ], p1 = 0.5, p0 = 0.5)
  expect_equal(as.character(screened$excluded_by_check), "3")
  expect_equal(as.character(screened$region_wise_remaining_regions), "2")
  expect_equal(screened$region_wise_estimate, retained$estimate, tolerance = 1e-8)
  expect_equal(screened$region_wise_se, retained$se_eif, tolerance = 1e-8)
  expect_equal(screened$region_wise_ci, retained$ci, tolerance = 1e-8)
})

test_that("omitting standard errors preserves the ATE point estimate", {
  dat <- make_ate_data()
  a <- regATE_dr(dat, p1 = 0.5, p0 = 0.5)
  b <- regATE_dr(dat, p1 = 0.5, p0 = 0.5, compute_se = FALSE)
  expect_equal(b$estimate, a$estimate)
  expect_null(b$se_eif)
  expect_null(b$ci)
})
