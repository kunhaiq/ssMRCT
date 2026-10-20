test_that("random-effects counts attain power with the smallest control count", {
  f <- c(0.2, 0.3, 0.5)
  sd_t <- c(1, 1.2, 0.9)
  sd_c <- c(1.1, 1, 1.2)
  ratio <- 2
  for (power in c(0.8, 0.9)) {
    ans <- ssCPnorm(0.025, 1 - power, 0.3, 0.08, sd_t, sd_c,
                    randRatio = ratio, target = 2, f = f)
    required <- (qnorm(0.975) + qnorm(power))^2 / 0.3^2
    information <- function(n) {
      sum(1 / (0.08^2 + (sd_t^2 / ratio + sd_c^2) / (n * f)))
    }
    expect_gte(information(ans$n0), required)
    expect_lt(information(ans$n0 - 1), required)
    expect_equal(ans$n1, ceiling(ratio * ans$n0))
    expect_equal(ans$n, ans$n0 + ans$n1)
  }
})

test_that("relabeling regions preserves the selected-region random-effects result", {
  f <- c(0.2, 0.3, 0.5)
  st <- c(1, 1.2, 0.9)
  sc <- c(1.1, 1, 1.2)
  order <- c(3, 1, 2)
  a <- ssCPnorm(0.025, 0.2, 0.3, 0.08, st, sc, target = 2, f = f)
  b <- ssCPnorm(0.025, 0.2, 0.3, 0.08, st[order], sc[order],
                target = match(2, order), f = f[order])
  expect_equal(a, b, tolerance = 1e-8)
})
