test_that("regional fraction solvers attain common CP targets", {
  for (target in c(0.8, 0.9)) {
    one <- regFrac(0.025, 0.9, CP = target, d = 0.5, sigmaTrt = 2)
    expect_equal(conProb(0.025, 0.9, rF = one$rF, d = 0.5,
                         sigmaTrt = 2)$CP, target, tolerance = 1e-6)
    two <- regFrac2(0.025, 0.9, CP = target)
    expect_equal(conProb2(0.025, 0.9, rF1 = two, d1 = 0.5,
                          sigmaTrt1 = 2)$CP, target, tolerance = 1e-6)
    expect_lt(two, one$rF)
  }
})

test_that("swapping independent trials preserves pooled Method I CP", {
  a <- conProb2(0.025, 0.8, 0.9, rF1 = 0.2, rF2 = 0.35,
                d1 = 1, d2 = 0.8, muCtrl1 = 0, muCtrl2 = 0.5,
                sigmaTrt1 = 4, sigmaCtrl1 = 3.5,
                sigmaTrt2 = 3.5, sigmaCtrl2 = 3,
                randRatio1 = 1, randRatio2 = 2)
  b <- conProb2(0.025, 0.9, 0.8, rF1 = 0.35, rF2 = 0.2,
                d1 = 0.8, d2 = 1, muCtrl1 = 0.5, muCtrl2 = 0,
                sigmaTrt1 = 3.5, sigmaCtrl1 = 3,
                sigmaTrt2 = 4, sigmaCtrl2 = 3.5,
                randRatio1 = 2, randRatio2 = 1)
  expect_equal(a$CP, b$CP, tolerance = 1e-6)
  expect_equal(a$N1, b$N2)
  expect_equal(a$N2, b$N1)
})

test_that("regional ordering does not change all-region Method II CP", {
  f1 <- c(0.2, 0.3, 0.5)
  f2 <- c(0.3, 0.4, 0.3)
  order <- c(3, 1, 2)
  set.seed(110)
  a <- conProbII2(0.025, 0.8, power2 = 0.9, rF1 = f1, rF2 = f2,
                  d1 = 1, d2 = 0.8, sigmaTrt1 = 4, randRatio2 = 2)
  set.seed(110)
  b <- conProbII2(0.025, 0.8, power2 = 0.9,
                  rF1 = f1[order], rF2 = f2[order],
                  d1 = 1, d2 = 0.8, sigmaTrt1 = 4, randRatio2 = 2)
  # The joint normal integral is numerical, so do not require exact equality.
  expect_equal(a$CP, b$CP, tolerance = 1e-4)
  expect_equal(a[c("N1", "N2")], b[c("N1", "N2")])
})
