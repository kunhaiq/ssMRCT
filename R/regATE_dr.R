#' Compute the doubly robust (DR) estimate and EIF-based standard error
#'
#' @description
#' Internal helper that implements the doubly robust estimator for the regional
#' average treatment effect \eqn{\tau} in the target region, together with the
#' efficient influence function (EIF) based standard error. The function fits
#' outcome, region-propensity and treatment-propensity models, optionally
#' adjusts for heterogeneous variance through the variance ratios \code{v1} and
#' \code{v0}, and returns all nuisance components for downstream reuse.
#'
#' @param dat A data frame with columns \code{y}, \code{trt}, \code{D}
#'   (target-region indicator) and covariates.
#' @param q A scalar giving the target-region proportion \eqn{\mathrm{Pr}(D=1)}.
#' @param covar_cols A character vector of covariate names.
#' @param compute_se Logical. Whether to compute the EIF-based standard error.
#' @param p1 Numeric scalar, numeric vector, or \code{NULL}. A constant treatment
#'   propensity \eqn{p_1(X)}; if \code{NULL} it is estimated by logistic regression
#'   within the target region.
#' @param p0 Same as \code{p1}, but for \eqn{p_0(X)}.
#' @param pi Numeric scalar, numeric vector, or \code{NULL}. A constant region
#'   propensity; if \code{NULL} it is estimated by logistic regression.
#' @param v1 Controls \eqn{v_1(X) = \mathrm{var}(Y \mid D=1, T=1) /
#'   \mathrm{var}(Y \mid D=0, T=1)}:
#'   \itemize{
#'     \item a numeric scalar: a fixed constant for all observations;
#'     \item \code{"const"}: a single constant estimated from the data;
#'     \item \code{"log_linear"}: a log-linear model is fitted for the conditional
#'           variance and \eqn{v_1(X)} is predicted pointwise.
#'   }
#' @param v0 Same as \code{v1}, but for \eqn{v_0(X) = \mathrm{var}(Y \mid D=1, T=0) /
#'   \mathrm{var}(Y \mid D=0, T=0)}.
#'
#' @return A list with components:
#' \item{tau_hat}{The doubly robust point estimate.}
#' \item{se}{The EIF-based standard error (if \code{compute_se = TRUE}).}
#' \item{m1_model, m0_model}{Fitted outcome models on the treated/control.}
#' \item{m1_pred, m0_pred}{Predicted outcome on all observations.}
#' \item{p1_model, p0_model, pi_model}{Fitted propensity models (or \code{NULL}).}
#' \item{p1_pred, p0_pred, pi_pred}{Predicted propensity scores.}
#' \item{v1, v0}{Variance ratio vectors.}
#' \item{v1_model, v0_model}{Fitted log-linear variance models (or \code{NULL}).}
#'
#' @keywords internal
#' @noRd
compute_dr_estimate <- function(dat, q, covar_cols, compute_se = TRUE,
                                p1 = NULL, p0 = NULL, pi = NULL,
                                v1 = "const", v0 = "const") {

  eif_se <- function(IF_vec) {
    sqrt(mean((IF_vec - mean(IF_vec))^2) / length(IF_vec))
  }

  n_all        <- nrow(dat)
  idx_D1       <- which(dat$D == 1)
  idx_D0       <- which(dat$D == 0)
  idx_trt_all  <- which(dat$trt == 1)
  idx_ctrl_all <- which(dat$trt == 0)

  covar_formula <- paste(covar_cols, collapse = " + ")

  # ---- Outcome models ----
  m1_formula <- as.formula(paste("y ~", covar_formula))
  m0_formula <- as.formula(paste("y ~", covar_formula))

  m1_full <- lm(m1_formula, data = dat[idx_trt_all, , drop = FALSE])
  m0_full <- lm(m0_formula, data = dat[idx_ctrl_all, , drop = FALSE])
  m1_full_pred <- predict(m1_full, newdata = dat)
  m0_full_pred <- predict(m0_full, newdata = dat)
  delta_full    <- m1_full_pred - m0_full_pred
  R1_full       <- dat$y - m1_full_pred
  R0_full       <- dat$y - m0_full_pred

  # ---- Region propensity pi(X) ----
  if (is.null(pi)) {
    region_ps_formula <- as.formula(paste("D ~", covar_formula))
    region_ps_model <- suppressWarnings(
      glm(region_ps_formula, data = dat, family = binomial)
    )
    pi_D1 <- predict(region_ps_model, newdata = dat, type = "response")
  } else if (length(pi) == 1) {
    pi_D1 <- rep(pi, n_all)
  } else {
    pi_D1 <- pi
  }

  # ---- Treatment propensity models p_d(X) ----
  dat_D1 <- dat[idx_D1, , drop = FALSE]
  dat_D0 <- dat[idx_D0, , drop = FALSE]

  # p1(X)
  if (is.null(p1)) {
    trt_formula <- as.formula(paste("trt ~", covar_formula))
    trt_ps_D1_model <- tryCatch(
      glm(trt_formula, data = dat_D1, family = binomial),
      error = function(e) NULL
    )
    if (!is.null(trt_ps_D1_model)) {
      pi_T_D1 <- predict(trt_ps_D1_model, newdata = dat, type = "response")
    } else {
      pi_T_D1 <- rep(0.5, n_all)
    }
  } else if (length(p1) == 1) {
    pi_T_D1 <- rep(p1, n_all)
  } else {
    pi_T_D1 <- p1
  }

  # p0(X)
  if (is.null(p0)) {
    trt_formula <- as.formula(paste("trt ~", covar_formula))
    trt_ps_D0_model <- tryCatch(
      glm(trt_formula, data = dat_D0, family = binomial),
      error = function(e) NULL
    )
    if (!is.null(trt_ps_D0_model)) {
      pi_T_D0 <- predict(trt_ps_D0_model, newdata = dat, type = "response")
    } else {
      pi_T_D0 <- rep(0.5, n_all)
    }
  } else if (length(p0) == 1) {
    pi_T_D0 <- rep(p0, n_all)
  } else {
    pi_T_D0 <- p0
  }

  # ---- Variance ratios v1, v0 ----
  # Helper: estimate a single constant ratio
  est_const_v1 <- function() {
    var_Y1_D1 <- var(dat$y[dat$D == 1 & dat$trt == 1])
    var_Y1_D0 <- var(dat$y[dat$D == 0 & dat$trt == 1])
    if (length(dat$y[dat$D == 1 & dat$trt == 1]) > 1 &&
        length(dat$y[dat$D == 0 & dat$trt == 1]) > 1 &&
        var_Y1_D0 > 0) var_Y1_D1 / var_Y1_D0 else 1
  }
  est_const_v0 <- function() {
    var_Y0_D1 <- var(dat$y[dat$D == 1 & dat$trt == 0])
    var_Y0_D0 <- var(dat$y[dat$D == 0 & dat$trt == 0])
    if (length(dat$y[dat$D == 1 & dat$trt == 0]) > 1 &&
        length(dat$y[dat$D == 0 & dat$trt == 0]) > 1 &&
        var_Y0_D0 > 0) var_Y0_D1 / var_Y0_D0 else 1
  }

  # v1
  v1_model <- NULL
  if (is.numeric(v1) && length(v1) == 1) {
    # Mode 1: fixed constant
    v1_vec <- rep(v1, n_all)
  } else if (is.character(v1) && v1 == "const") {
    # Mode 2: estimate constant
    v1_vec <- rep(est_const_v1(), n_all)
  } else if (is.character(v1) && v1 == "log_linear") {
    # Mode 3: log-linear model with D interaction (on T=1 only)
    idx_T1 <- which(dat$trt == 1)
    interact_terms <- paste(paste0(covar_cols, ":D"), collapse = " + ")
    logvar_formula <- as.formula(paste(
      "log(pmax(resid_sq, 1e-8)) ~", covar_formula, "+ D +", interact_terms))
    if (length(idx_T1) > 2 * length(covar_cols) + 3 &&
        length(unique(dat$D[idx_T1])) == 2) {
      df_T1 <- dat[idx_T1, , drop = FALSE]
      df_T1$resid_sq <- (df_T1$y - m1_full_pred[idx_T1])^2
      fit_logvar_T1 <- lm(logvar_formula, data = df_T1)
      dat_D1 <- dat; dat_D1$D <- 1
      dat_D0 <- dat; dat_D0$D <- 0
      pred_D1 <- predict(fit_logvar_T1, newdata = dat_D1)
      pred_D0 <- predict(fit_logvar_T1, newdata = dat_D0)
      v1_vec <- exp(pred_D1 - pred_D0)
      v1_model <- fit_logvar_T1
    } else {
      v1_vec <- rep(est_const_v1(), n_all)
    }
  } else {
    stop("v1 must be a numeric scalar, 'const', or 'log_linear'")
  }

  # v0
  v0_model <- NULL
  if (is.numeric(v0) && length(v0) == 1) {
    v0_vec <- rep(v0, n_all)
  } else if (is.character(v0) && v0 == "const") {
    v0_vec <- rep(est_const_v0(), n_all)
  } else if (is.character(v0) && v0 == "log_linear") {
    # Mode 3: log-linear model with D interaction (on T=0 only)
    idx_T0 <- which(dat$trt == 0)
    interact_terms <- paste(paste0(covar_cols, ":D"), collapse = " + ")
    logvar_formula <- as.formula(paste(
      "log(pmax(resid_sq, 1e-8)) ~", covar_formula, "+ D +", interact_terms))
    if (length(idx_T0) > 2 * length(covar_cols) + 3 &&
        length(unique(dat$D[idx_T0])) == 2) {
      df_T0 <- dat[idx_T0, , drop = FALSE]
      df_T0$resid_sq <- (df_T0$y - m0_full_pred[idx_T0])^2
      fit_logvar_T0 <- lm(logvar_formula, data = df_T0)
      dat_D1 <- dat; dat_D1$D <- 1
      dat_D0 <- dat; dat_D0$D <- 0
      pred_D1 <- predict(fit_logvar_T0, newdata = dat_D1)
      pred_D0 <- predict(fit_logvar_T0, newdata = dat_D0)
      v0_vec <- exp(pred_D1 - pred_D0)
      v0_model <- fit_logvar_T0
    } else {
      v0_vec <- rep(est_const_v0(), n_all)
    }
  } else {
    stop("v0 must be a numeric scalar, 'const', or 'log_linear'")
  }

  # ---- Weights ----
  W1_num   <- dat$D * dat$trt * pi_D1 +
              (1 - dat$D) * dat$trt * pi_D1 * v1_vec
  W1_denom <- pi_D1 * pi_T_D1 + (1 - pi_D1) * pi_T_D0 * v1_vec
  W1 <- W1_num / W1_denom

  W0_num   <- dat$D * (1 - dat$trt) * pi_D1 +
              (1 - dat$D) * (1 - dat$trt) * pi_D1 * v0_vec
  W0_denom <- pi_D1 * (1 - pi_T_D1) + (1 - pi_D1) * (1 - pi_T_D0) * v0_vec
  W0 <- W0_num / W0_denom

  # ---- DR estimate ----
  tau_hat <- (1 / q) * mean(dat$D * delta_full + W1 * R1_full - W0 * R0_full)

  # ---- EIF-based SE ----
  se_val <- NULL
  if (compute_se) {
    IF_vec <- (1 / q) * (dat$D * (delta_full - tau_hat) +
                         W1 * R1_full - W0 * R0_full)
    se_val <- eif_se(IF_vec)
  }

  p1_is_const <- !is.null(p1) && length(p1) == 1
  p0_is_const <- !is.null(p0) && length(p0) == 1
  pi_is_const <- !is.null(pi) && length(pi) == 1

  list(
    tau_hat     = tau_hat,
    se          = se_val,
    m1_model    = m1_full,
    m0_model    = m0_full,
    m1_pred     = m1_full_pred,
    m0_pred     = m0_full_pred,
    p1_model    = if (p1_is_const) NULL else trt_ps_D1_model,
    p1_pred     = pi_T_D1,
    p0_model    = if (p0_is_const) NULL else trt_ps_D0_model,
    p0_pred     = pi_T_D0,
    pi_model    = if (pi_is_const) NULL else region_ps_model,
    pi_pred     = pi_D1,
    v1          = v1_vec,
    v0          = v0_vec,
    v1_model    = v1_model,
    v0_model    = v0_model
  )
}

#' @name regATE_dr
#'
#' @title Doubly robust estimator for regional average treatment effects
#'
#' @description
#' This function implements a doubly robust estimator for the average treatment
#' effect (ATE) in a target region within a Multi-Regional Clinical Trial (MRCT).
#' It optionally performs region-wise ANCOVA F-test screening to adaptively
#' borrow information from source regions that are compatible with the target
#' region.
#'
#' @param dat A data frame with columns \code{y}, \code{trt}, \code{region} and
#'   numeric covariates.
#' @param target_region Integer. The region of interest (default \code{1}).
#' @param compute_se Logical. Whether to compute the EIF-based standard error
#'   (default \code{TRUE}).
#' @param region_wise Logical. Whether to run the region-wise ANCOVA F-test
#'   screening (default \code{FALSE}).
#' @param p1 Numeric scalar, numeric vector, or \code{NULL}. A constant treatment
#'   propensity \eqn{p_1(X)}; if \code{NULL} it is estimated by logistic
#'   regression within the target region.
#' @param p0 Same as \code{p1}, but for \eqn{p_0(X)}.
#' @param pi Numeric scalar, numeric vector, or \code{NULL}. A constant region
#'   propensity; if \code{NULL} it is estimated by logistic regression.
#' @param v1 Controls the variance ratio \eqn{v_1(X)}: a numeric scalar (fixed
#'   constant), \code{"const"} (default; estimated constant), or
#'   \code{"log_linear"} (pointwise prediction from a log-linear variance model).
#' @param v0 Same as \code{v1}, but for \eqn{v_0(X)}.
#' @param config A named list to override screening settings. The only supported
#'   element is \code{check_p_threshold} (default \code{0.05}), the p-value
#'   threshold used in the region-wise ANCOVA F-test.
#'
#' @details
#' The doubly robust estimator is given by
#' \deqn{
#' \widehat{\tau}_{\mathrm{dr}} = \frac{1}{q} \sum_{i=1}^{n} \left[
#'   D_i\{\widehat{m}_1(X_i) - \widehat{m}_0(X_i)\} +
#'   \widehat{W}_{1,i}\{Y_i - \widehat{m}_1(X_i)\} -
#'   \widehat{W}_{0,i}\{Y_i - \widehat{m}_0(X_i)\}
#' \right],
#' }
#' where \eqn{q = \mathrm{Pr}(D=1)} is the target-region proportion and
#' \eqn{\widehat{W}_{1,i}}, \eqn{\widehat{W}_{0,i}} are weights constructed from
#' the region and treatment propensity scores (optionally adjusted by the
#' variance ratios \eqn{v_1(X), v_0(X)}). When \code{region_wise = TRUE}, each
#' source region is screened by an ANCOVA F-test comparing the reduced model
#' \eqn{y \sim X} with the full model \eqn{y \sim X + D_s + X:D_s}; regions that
#' fail the test (p-value below \code{check_p_threshold} in either the treated or
#' the control group) are excluded from the borrowing pool.
#'
#' @returns A list of class \code{regATE_dr} with components:
#' \describe{
#'   \item{\code{estimate}}{The simple-pooling doubly robust estimate.}
#'   \item{\code{se_eif}}{The EIF-based standard error (if \code{compute_se = TRUE}).}
#'   \item{\code{ci}}{The 95% Wald confidence interval (if \code{compute_se = TRUE}).}
#'   \item{\code{remaining_regions}}{The candidate source regions.}
#'   \item{\code{config}}{The (possibly merged) configuration list.}
#'   \item{\code{check_results}}{A data frame of region-wise ANCOVA F-test
#'         results (only if \code{region_wise = TRUE}).}
#'   \item{\code{excluded_by_check}}{Source regions excluded by the screening
#'         (only if \code{region_wise = TRUE}).}
#'   \item{\code{region_wise_estimate}, \code{region_wise_se},
#'         \code{region_wise_ci}}{The estimate, standard error and confidence
#'         interval based on the retained regions (only if
#'         \code{region_wise = TRUE}).}
#'   \item{\code{region_wise_remaining_regions}}{Source regions retained after
#'         screening (only if \code{region_wise = TRUE}).}
#' }
#'
#' @examples
#' set.seed(42)
#' n <- 600
#' region <- c(rep(1, 150), rep(2, 200), rep(3, 250))
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' trt <- rbinom(n, 1, 0.5)
#' y <- x1 + x2 + trt * (1 + 0.5 * x1) + rnorm(n)
#' dat <- data.frame(y = y, trt = trt, region = region, x1 = x1, x2 = x2)
#'
#' ## Simple-pooling doubly robust estimate
#' res <- regATE_dr(dat, target_region = 1)
#' print(res)
#'
#' ## Region-wise screening
#' res2 <- regATE_dr(dat, target_region = 1, region_wise = TRUE)
#' print(res2)
#'
#' @importFrom stats as.formula binomial glm lm predict qnorm var anova
#' @importFrom utils modifyList
#' @export
regATE_dr <- function(
    dat,
    target_region = 1,
    compute_se    = TRUE,
    region_wise   = FALSE,
    p1            = NULL,
    p0            = NULL,
    pi            = NULL,
    v1            = "const",
    v0            = "const",
    config        = list()
) {

  K <- length(unique(dat$region))

  default_config <- list(
    check_p_threshold = 0.05
  )
  config <- modifyList(default_config, config)

  required_cols <- c("y", "trt", "region")
  missing_cols  <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("dat is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!(target_region %in% unique(dat$region))) {
    stop("target_region = ", target_region, " not found in dat$region")
  }

  if (length(unique(dat$region)) < 2) {
    stop("dat must contain at least 2 regions")
  }

  all_cols   <- names(dat)
  covar_cols <- setdiff(all_cols, c("y", "trt", "region"))
  covar_cols <- covar_cols[sapply(dat[covar_cols], is.numeric)]

  if (length(covar_cols) == 0) {
    stop("No numeric covariates found in dat")
  }

  all_regions       <- sort(unique(dat$region))
  source_candidates <- setdiff(all_regions, target_region)

  dat$D <- as.integer(dat$region == target_region)
  q     <- mean(dat$D)

  # ---- Simple pooling estimate ----
  idx_pool     <- which(dat$region %in% c(target_region, source_candidates))
  dat_pool     <- dat[idx_pool, , drop = FALSE]
  dat_pool$D   <- as.integer(dat_pool$region == target_region)
  q_pool       <- mean(dat_pool$D)

  pool_est <- compute_dr_estimate(dat_pool, q_pool, covar_cols, compute_se = compute_se,
                                  p1 = p1, p0 = p0, pi = pi, v1 = v1, v0 = v0)

  out <- list(
    estimate              = pool_est$tau_hat,
    se_eif                = pool_est$se,
    ci                    = if (compute_se && !is.null(pool_est$se)) {
      c(pool_est$tau_hat - qnorm(0.975) * pool_est$se,
        pool_est$tau_hat + qnorm(0.975) * pool_est$se)
    } else NULL,
    remaining_regions = source_candidates,
    config                = config
  )

  if (!region_wise) {
    class(out) <- "regATE_dr"
    return(out)
  }

  # ---- Step 2: Region-wise checks ----
  excluded_by_check <- character(0)
  check_rows        <- list()

  for (src in source_candidates) {
    idx_src   <- which(dat$region %in% c(target_region, src))
    dat_check <- dat[idx_src, , drop = FALSE]
    dat_check$D_src <- as.integer(dat_check$region == src)

    dat_trt <- dat_check[dat_check$trt == 1, , drop = FALSE]
    p_trt   <- NA

    if (nrow(dat_trt) > length(covar_cols) + 2 &&
        length(unique(dat_trt$D_src)) == 2) {

      covar_formula <- paste(covar_cols, collapse = " + ")
      reduced_formula <- as.formula(paste("y ~", covar_formula))
      full_formula <- as.formula(
        paste("y ~", covar_formula, "+ D_src +",
              paste(paste0(covar_cols, ":D_src"), collapse = " + "))
      )

      fit_reduced <- tryCatch(lm(reduced_formula, data = dat_trt), error = function(e) NULL)
      fit_full    <- tryCatch(lm(full_formula,    data = dat_trt), error = function(e) NULL)

      if (!is.null(fit_reduced) && !is.null(fit_full)) {
        aov_result <- anova(fit_reduced, fit_full)
        p_trt <- aov_result[2, "Pr(>F)"]
      }
    }

    dat_ctrl <- dat_check[dat_check$trt == 0, , drop = FALSE]
    p_ctrl   <- NA

    if (nrow(dat_ctrl) > length(covar_cols) + 2 &&
        length(unique(dat_ctrl$D_src)) == 2) {

      covar_formula <- paste(covar_cols, collapse = " + ")
      reduced_formula <- as.formula(paste("y ~", covar_formula))
      full_formula <- as.formula(
        paste("y ~", covar_formula, "+ D_src +",
              paste(paste0(covar_cols, ":D_src"), collapse = " + "))
      )

      fit_reduced <- tryCatch(lm(reduced_formula, data = dat_ctrl), error = function(e) NULL)
      fit_full    <- tryCatch(lm(full_formula,    data = dat_ctrl), error = function(e) NULL)

      if (!is.null(fit_reduced) && !is.null(fit_full)) {
        aov_result <- anova(fit_reduced, fit_full)
        p_ctrl <- aov_result[2, "Pr(>F)"]
      }
    }

    pass <- TRUE
    if (!is.na(p_trt) && p_trt < config$check_p_threshold) pass <- FALSE
    if (!is.na(p_ctrl) && p_ctrl < config$check_p_threshold) pass <- FALSE

    if (!pass) excluded_by_check <- c(excluded_by_check, as.character(src))

    check_rows[[length(check_rows) + 1]] <- data.frame(
      source_region   = src,
      p_value_trt     = p_trt,
      p_value_ctrl    = p_ctrl,
      passed          = pass,
      stringsAsFactors = FALSE
    )
  }

  check_results <- do.call(rbind, check_rows)
  source_pool   <- setdiff(source_candidates, excluded_by_check)

  out$check_results     <- check_results
  out$excluded_by_check <- excluded_by_check

  if (length(source_pool) > 0) {
    idx_wf     <- which(dat$region %in% c(target_region, source_pool))
    dat_wf     <- dat[idx_wf, , drop = FALSE]
    dat_wf$D   <- as.integer(dat_wf$region == target_region)
    q_wf       <- mean(dat_wf$D)
    wf_est <- compute_dr_estimate(dat_wf, q_wf, covar_cols, compute_se = compute_se,
                                  p1 = p1, p0 = p0, pi = pi, v1 = v1, v0 = v0)
  } else {
    idx_wf     <- which(dat$region == target_region)
    dat_wf     <- dat[idx_wf, , drop = FALSE]
    dat_wf$D   <- 1
    q_wf       <- 1
    wf_est <- compute_dr_estimate(dat_wf, q_wf, covar_cols, compute_se = compute_se,
                                  p1 = p1, p0 = p0, pi = pi, v1 = v1, v0 = v0)
  }

  out$region_wise_estimate      <- wf_est$tau_hat
  out$region_wise_se            <- wf_est$se
  out$region_wise_ci            <- if (compute_se && !is.null(wf_est$se)) {
    c(wf_est$tau_hat - qnorm(0.975) * wf_est$se,
      wf_est$tau_hat + qnorm(0.975) * wf_est$se)
  } else NULL
  out$region_wise_remaining_regions <- source_pool

  class(out) <- "regATE_dr"
  return(out)
}

#' @rdname regATE_dr
#' @param x An object of class \code{regATE_dr}.
#' @param ... Additional arguments passed to other methods (currently unused).
#' @export
print.regATE_dr <- function(x, ...) {
  cat("\n")
  cat("============================================================\n")
  cat("  regATE_dr: Doubly Robust Regional ATE Estimation\n")
  cat("============================================================\n\n")

  cat("-- Simple pooling estimate (all candidate sources) --\n")
  cat(sprintf("  Remaining regions     : %s\n",
              paste(x$remaining_regions, collapse = ", ")))
  cat(sprintf("  Estimate              : %.4f\n", x$estimate))
  if (!is.null(x$se_eif)) {
    cat(sprintf("  SE (EIF)              : %.4f\n", x$se_eif))
    cat(sprintf("  95%% CI                : [%.4f, %.4f]\n", x$ci[1], x$ci[2]))
  }

  has_region_wise <- !is.null(x$region_wise_estimate) || !is.null(x$check_results)

  if (has_region_wise) {
    cat("\n")
    cat("============================================================\n")
    cat("  Workflow-based estimation\n")
    cat("============================================================\n\n")

    cat("--- Region-wise ANCOVA F-test ---\n")
    cat(sprintf("  Threshold: p < %.2f\n", x$config$check_p_threshold))
    cat("  Full model   : y ~ X + D_src + X:D_src\n")
    cat("  Reduced model: y ~ X\n")
    if (!is.null(x$check_results)) {
      print(x$check_results, row.names = FALSE)
    }
    if (length(x$excluded_by_check) > 0) {
      cat(sprintf("  => Excluded by check: %s\n",
                  paste(x$excluded_by_check, collapse = ", ")))
    } else {
      cat("  => All candidate source regions passed.\n")
    }
    if (!is.null(x$region_wise_estimate)) {
      cat(sprintf("  Region-wise estimate   : %.4f\n", x$region_wise_estimate))
      cat(sprintf("  Remaining regions      : %s\n",
                  if (length(x$region_wise_remaining_regions) > 0)
                    paste(x$region_wise_remaining_regions, collapse = ", ")
                  else "(none - target region only)"))
      if (!is.null(x$region_wise_se)) {
        cat(sprintf("  Region-wise SE         : %.4f\n", x$region_wise_se))
        cat(sprintf("  Region-wise 95%% CI     : [%.4f, %.4f]\n",
                    x$region_wise_ci[1], x$region_wise_ci[2]))
      }
    }
  }

  cat("\n")
  invisible(x)
}