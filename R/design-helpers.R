# Shared validation and design calculations (internal).
.check_number <- function(x, name, lower = -Inf, upper = Inf,
                          lower_open = FALSE, upper_open = FALSE, scalar = TRUE) {
  if (!is.numeric(x) || !length(x) || (scalar && length(x) != 1L) ||
      any(!is.finite(x)) || any(if (lower_open) x <= lower else x < lower) ||
      any(if (upper_open) x >= upper else x > upper)) {
    stop(name, " is outside its valid numeric range", call. = FALSE)
  }
  invisible(x)
}
.check_design <- function(alpha, power, pi = 0.5) {
  .check_number(alpha, "alpha", 0, 0.5, TRUE, TRUE)
  .check_number(power, "power", alpha, 1, TRUE, TRUE)
  .check_number(pi, "pi", 0, 1, FALSE, TRUE)
}
.fractions <- function(f, name = "f") {
  .check_number(f, name, 0, 1, TRUE, FALSE, FALSE)
  if (abs(sum(f) - 1) > 1e-10) stop(name, " must sum to 1", call. = FALSE)
  f / sum(f)
}
.regional <- function(x, k, name, lower = 0, upper = Inf, open = TRUE) {
  .check_number(x, name, lower, upper, open, FALSE, FALSE)
  if (!length(x) %in% c(1L, k)) stop(name, " must be scalar or have one value per region", call. = FALSE)
  rep(x, length.out = k)
}
.fixed_design <- function(alpha, power, d, sigmaTrt, sigmaCtrl, ratio) {
  .check_design(alpha, power)
  .check_number(d, "d", 0, lower_open = TRUE)
  .check_number(sigmaTrt, "sigmaTrt", 0, lower_open = TRUE)
  .check_number(sigmaCtrl, "sigmaCtrl", 0, lower_open = TRUE)
  .check_number(ratio, "randRatio", 0, lower_open = TRUE)
  nc <- (sigmaTrt^2 / ratio + sigmaCtrl^2) * (qnorm(1-alpha) + qnorm(power))^2 / d^2
  .check_number(nc, "continuous control sample size", 0, lower_open = TRUE)
  list(nc = nc, nt = ratio * nc, nCtrl = ceiling(nc), nTrt = ceiling(ratio * ceiling(nc)))
}
# A degenerate normal variable is nonnegative iff its mean is nonnegative.
.normal_nonnegative <- function(mu, sd) {
  if (sd == 0) as.numeric(mu >= 0) else pnorm(mu / sd)
}
.method1_cp <- function(alpha, power, pi, fraction, studies = 1L) {
  if (fraction == 0) return(0.5) # limiting value; not an executable allocation
  if (fraction == 1) return(1)
  z <- qnorm(1-alpha) + qnorm(power)
  den <- sqrt(studies * (1/fraction - 1))
  if (studies == 1L) {
    return(integrate(function(u) pnorm((1-pi)*(u+z)/den)*dnorm(u),
                     -qnorm(power), Inf)$value / power)
  }
  integrate(function(v) vapply(v, function(vi) {
    integrate(function(u) pnorm((1-pi)*(u+vi+2*z)/den)*dnorm(u),
              -qnorm(power), Inf)$value * dnorm(vi)
  }, numeric(1)), -qnorm(power), Inf)$value / power^2
}
.binary_design <- function(alpha, power, pTrt, pCtrl, ratio, f) {
  .check_number(pTrt, "pTrt", 0, 1, TRUE, TRUE)
  .check_number(pCtrl, "pCtrl", 0, 1, TRUE, TRUE)
  x <- .fixed_design(alpha, power, pTrt-pCtrl, sqrt(pTrt*(1-pTrt)),
                     sqrt(pCtrl*(1-pCtrl)), ratio)
  list(nTrtReg = ceiling(f*x$nTrt), nCtrlReg = ceiling(f*x$nCtrl),
       pTrt = pTrt, pCtrl = pCtrl)
}
# Simulate regional sufficient statistics in bounded batches, not individual outcomes.
.binary_cp <- function(designs, alpha, B) {
  .check_number(B, "B", 1, .Machine$integer.max)
  if (B != floor(B)) stop("B must be an integer", call. = FALSE)
  k <- length(designs[[1]]$nTrtReg)
  totalT <- Reduce(`+`, lapply(designs, `[[`, "nTrtReg"))
  totalC <- Reduce(`+`, lapply(designs, `[[`, "nCtrlReg"))
  successCount <- consistencyCount <- 0
  done <- 0
  while (done < B) {
    b <- min(10000, B-done)
    pooledT <- pooledC <- matrix(0, b, k)
    success <- rep(TRUE, b)
    for (des in designs) {
      yt <- vapply(des$nTrtReg, function(n) rbinom(b, n, des$pTrt), numeric(b))
      yc <- vapply(des$nCtrlReg, function(n) rbinom(b, n, des$pCtrl), numeric(b))
      dim(yt) <- dim(yc) <- c(b, k)
      pt <- rowSums(yt)/sum(des$nTrtReg)
      pc <- rowSums(yc)/sum(des$nCtrlReg)
      se <- sqrt(pt*(1-pt)/sum(des$nTrtReg) + pc*(1-pc)/sum(des$nCtrlReg))
      # This also defines zero-variance cases without 0/0: a zero effect never rejects.
      success <- success & (pt-pc > qnorm(1-alpha)*se)
      pooledT <- pooledT + yt
      pooledC <- pooledC + yc
    }
    consistent <- rowSums(sweep(pooledT, 2, totalC, `*`) >=
                           sweep(pooledC, 2, totalT, `*`)) == k
    successCount <- successCount + sum(success)
    consistencyCount <- consistencyCount + sum(success & consistent)
    done <- done + b
  }
  cp <- if (successCount) consistencyCount/successCount else NA_real_
  if (!successCount) warning("No overall successes: conditional CP cannot be estimated", call. = FALSE)
  list(CP = cp, achievedPower = successCount/B,
       MCSE = if (successCount) sqrt(cp*(1-cp)/successCount) else NA_real_,
       successes = successCount, consistentSuccesses = consistencyCount, B = B)
}
.random_design <- function(alpha, beta, delta, tau, Omega, pi, ratio, target, f) {
  .check_number(beta, "beta", 0, 1, TRUE, TRUE)
  .check_design(alpha, 1-beta, pi)
  .check_number(delta, "delta", 0, lower_open = TRUE)
  .check_number(tau, "tau", 0)
  .check_number(ratio, "randRatio", 0, lower_open = TRUE)
  f <- .fractions(f)
  if (length(f) < 2) stop("At least two regions are required", call.=FALSE)
  .check_number(target, "target", 1, length(f))
  if (target != floor(target)) stop("target must be an integer", call. = FALSE)
  Omega <- .regional(Omega, length(f), "Omega")
  z <- qnorm(1-alpha) + qnorm(1-beta)
  rhs <- z^2/delta^2
  if (tau^2*rhs >= length(f)) stop("No finite sample size attains the requested power", call. = FALSE)
  information <- function(n) sum(n*f/(tau^2*n*f + Omega))
  hi <- 1
  while (information(hi) < rhs) {
    hi <- hi*2
    if (!is.finite(hi) || hi > 2^52) stop("Required sample size exceeds numerical precision", call. = FALSE)
  }
  nContinuous <- uniroot(function(n) information(n)-rhs, c(0, hi), tol=1e-9)$root
  n0 <- max(1, ceiling(nContinuous))
  # Guard the minimal integer solution against root-finding round-off.
  if (information(n0) < rhs) n0 <- n0+1
  if (n0 > 1 && information(n0-1) >= rhs) n0 <- n0-1
  # CP uses the same continuous design equation as the published approximation.
  shrink <- tau^2*nContinuous*f/(tau^2*nContinuous*f + Omega)
  den <- sqrt(shrink[target]*sum(shrink[-target]))
  cp <- integrate(function(u) .normal_nonnegative((1-pi)*(u+z), den)*dnorm(u),
                  qnorm(beta), Inf)$value/(1-beta)
  n1 <- ceiling(n0*ratio)
  list(n0=n0, n1=n1, n=n0+n1, CP=cp)
}
