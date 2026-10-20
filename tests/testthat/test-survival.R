test_that("homogeneous PH design agrees with the event-count formula", {
  alpha <- 0.025
  power <- 0.8
  hr <- 0.7
  rate <- 0.05
  follow_up <- 24
  for (ratio in c(1, 2)) {
    # Schoenfeld event requirement, then convert events to participants.
    events <- (qnorm(1 - alpha) + qnorm(power))^2 *
      (1 + ratio)^2 / (ratio * log(hr)^2)
    event_probability <- ((1 - exp(-rate * follow_up)) +
      ratio * (1 - exp(-rate * hr * follow_up))) / (1 + ratio)
    control <- ceiling(events / event_probability / (1 + ratio))
    ans <- ssCPsurvPH(alpha, 1 - power, rate, hr, randRatio = ratio,
                      target = 1, f = c(0.2, 0.3, 0.5), L = follow_up)
    expect_equal(ans$n0, control)
    expect_equal(ans$n1, ceiling(ratio * control))
    expect_equal(ans$CP, 1, tolerance = 1e-8)
  }
})

test_that("RMST variance agrees with an independent censored exponential integral", {
  eta <- 10
  reference_variance <- function(rate) {
    integrate(function(t) {
      tail_area <- (exp(-rate * t) - exp(-rate * eta)) / rate
      censor_survival <- pmin(1, (24 - t) / 22)
      tail_area^2 * rate / (exp(-rate * t) * censor_survival)
    }, 0, eta, rel.tol = 1e-9)$value
  }
  ans <- ssCPsurvNPH(0.025, 0.2,
    lambda0a = 0.1, lambda0b = 0.1, lambda1a = 0.07, lambda1b = 0.07,
    breakpoint = 5, unimin = 2, unimax = 24, eta = eta,
    target = 1, f = c(0.2, 0.3, 0.5))
  expect_equal(ans$sigma2.0r, rep(reference_variance(0.1), 3), tolerance = 1e-6)
  expect_equal(ans$sigma2.1r, rep(reference_variance(0.07), 3), tolerance = 1e-6)
  reference_effect <- (1 - exp(-0.07 * eta)) / 0.07 -
    (1 - exp(-0.1 * eta)) / 0.1
  expect_equal(ans$Dr, rep(reference_effect, 3), tolerance = 1e-7)
})
