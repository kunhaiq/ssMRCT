# Internal DR estimator; dat$D marks the target region and q is its fraction.
compute_dr_estimate <- function(dat, q, covar_cols, compute_se = TRUE,
                                p1 = NULL, p0 = NULL, pi = NULL,
                                v1 = "const", v0 = "const") {

  eif_se <- function(IF_vec) {
    sqrt(mean((IF_vec - mean(IF_vec))^2) / length(IF_vec))
  }

  trt_ps_D1_model <- trt_ps_D0_model <- region_ps_model <- NULL
  n_all        <- nrow(dat)
  idx_D1       <- which(dat$D == 1)
  idx_D0       <- which(dat$D == 0)
  idx_trt_all  <- which(dat$trt == 1)
  idx_ctrl_all <- which(dat$trt == 0)

  covar_formula <- paste(covar_cols, collapse = " + ")

  # Outcome models
  m1_formula <- as.formula(paste("y ~", covar_formula))
  m0_formula <- as.formula(paste("y ~", covar_formula))

  m1_full <- lm(m1_formula, data = dat[idx_trt_all, , drop = FALSE])
  m0_full <- lm(m0_formula, data = dat[idx_ctrl_all, , drop = FALSE])
  m1_full_pred <- predict(m1_full, newdata = dat)
  m0_full_pred <- predict(m0_full, newdata = dat)
  delta_full    <- m1_full_pred - m0_full_pred
  R1_full       <- dat$y - m1_full_pred
  R0_full       <- dat$y - m0_full_pred

  # Region propensity pi(X)
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

  # Treatment propensity models p_d(X)
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

  # Variance ratios v1, v0
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

  if (any(!is.finite(c(m1_full_pred, m0_full_pred, pi_D1, pi_T_D1, pi_T_D0,
                       v1_vec, v0_vec))) || any(v1_vec <= 0 | v0_vec <= 0)) {
    stop("Nuisance predictions must be finite, with positive variance ratios")
  }
  # Residual weights combine regional overlap, treatment assignment and variance.

  W1_num   <- dat$D * dat$trt * pi_D1 +
              (1 - dat$D) * dat$trt * pi_D1 * v1_vec
  W1_denom <- pi_D1 * pi_T_D1 + (1 - pi_D1) * pi_T_D0 * v1_vec
  if (any(W1_denom <= 0)) stop("Treatment weights are undefined; check overlap")
  W1 <- W1_num / W1_denom

  W0_num   <- dat$D * (1 - dat$trt) * pi_D1 +
              (1 - dat$D) * (1 - dat$trt) * pi_D1 * v0_vec
  W0_denom <- pi_D1 * (1 - pi_T_D1) + (1 - pi_D1) * (1 - pi_T_D0) * v0_vec
  if (any(W0_denom <= 0)) stop("Control weights are undefined; check overlap")
  W0 <- W0_num / W0_denom

  # DR estimate
  tau_hat <- (1 / q) * mean(dat$D * delta_full + W1 * R1_full - W0 * R0_full)

  # EIF-based SE
  se_val <- NULL
  if (compute_se) {
    IF_vec <- (1 / q) * (dat$D * (delta_full - tau_hat) +
                         W1 * R1_full - W0 * R0_full)
    se_val <- eif_se(IF_vec)
  }


  list(
    tau_hat     = tau_hat,
    se          = se_val,
    m1_model    = m1_full,
    m0_model    = m0_full,
    m1_pred     = m1_full_pred,
    m0_pred     = m0_full_pred,
    p1_model    = trt_ps_D1_model,
    p1_pred     = pi_T_D1,
    p0_model    = trt_ps_D0_model,
    p0_pred     = pi_T_D0,
    pi_model    = region_ps_model,
    pi_pred     = pi_D1,
    v1          = v1_vec,
    v0          = v0_vec,
    v1_model    = v1_model,
    v0_model    = v0_model
  )
}

#' Doubly robust estimation of a regional treatment effect
#'
#' @description Estimate the average treatment effect (ATE) in a target
#' region using covariate information from other regions. Optionally screen
#' source regions before borrowing information.
#'
#' @param dat A data frame with numeric outcome `y`, binary treatment `trt`
#'   (1 = treatment), region label `region`, and at least one numeric covariate.
#'   All other numeric columns are covariates. Values must be complete and finite;
#'   every region must contain both arms. Use valid R covariate names; `D`,
#'   `D_src` and `resid_sq` are reserved for internal use.
#' @param target_region A single label identifying the target region in `dat`.
#' @param compute_se Calculate the influence-function standard error and 95% interval.
#' @param region_wise Screen each source region using separate treatment-arm
#'   and control-arm ANCOVA tests before refitting.
#' @param p1,p0 Treatment assignment probabilities in the target and source
#'   populations. Supply `NULL` to fit logistic models, a scalar constant,
#'   or a vector aligned with all rows of `dat`. Supplied values must be in `(0, 1)`.
#' @param pi Probability of membership in the target region, given covariates.
#'   Accepts `NULL`, a scalar or a row-aligned vector, as for `p1`.
#'   This is a selection probability, not the CP preservation fraction.
#' @param v1,v0 Target-to-source variance ratios in the treatment and control
#'   arms. Use a positive scalar, `"const"` for a marginal variance ratio,
#'   or `"log_linear"` for a model based on log squared residuals.
#' @param config A named list with optional `check_p_threshold` in `(0, 1)`
#'   (default 0.05). Other entries are not supported.
#'
#' @details
#' Let \eqn{D_i} indicate target-region membership and
#' \eqn{\widehat q=n^{-1}\sum_i D_i}. The estimator is
#' \deqn{\widehat\tau = \frac{1}{n\widehat q}\sum_{i=1}^n
#'   \left[D_i\widehat\Delta(X_i)+\widehat W_{1i}\widehat U_{1i}
#'   -\widehat W_{0i}\widehat U_{0i}\right].}
#' Here \eqn{\widehat\Delta=\widehat m_1-\widehat m_0},
#' \eqn{\widehat U_t=Y-\widehat m_t(X)}, and the residual weights combine
#' selection scores, treatment assignment and variance ratios.
#' Outcome models are linear; estimated probability models are logistic.
#'
#' Borrowing assumes conditional outcome-mean exchangeability across the
#' regions retained. Double robustness concerns correct specification of
#' either outcome models or both selection and treatment-assignment models;
#' it does not remove the exchangeability assumption.
#'
#' @section Source screening:
#' Each arm compares `y ~ X` with `y ~ X + D_src + X:D_src`.
#' A source is excluded if either p-value is below `check_p_threshold`.
#' Treatment probability vectors are subset with the retained rows. Selection
#' probabilities are re-estimated after exclusions, because the population
#' changes. If no sources remain, the estimator uses only the target region.
#' The reported influence-function interval does not account for uncertainty
#' from selecting the source regions.
#'
#' @returns A list of class `regATE_dr` with:
#' \describe{
#'   \item{estimate}{ATE using all candidate source regions.}
#'   \item{se_eif, ci}{Standard error and 95% Wald interval, or `NULL` if disabled.}
#'   \item{remaining_regions}{All candidate source-region labels.}
#'   \item{config}{Screening settings used.}
#' }
#' With `region_wise = TRUE`, the list also contains:
#' \describe{
#'   \item{check_results}{Source labels, arm-specific p-values and pass indicators.}
#'   \item{excluded_by_check}{Excluded source labels.}
#'   \item{region_wise_estimate}{ATE after source screening.}
#'   \item{region_wise_se, region_wise_ci}{Its standard error and 95% interval.}
#'   \item{region_wise_remaining_regions}{Retained source labels.}
#' }
#' @examples
#' set.seed(42)
#' dat <- data.frame(region = rep(1:3, each = 100),
#'                   trt = rep(c(0, 1), 150), x = rnorm(300))
#' dat$y <- dat$x + dat$trt + rnorm(300)
#' regATE_dr(dat, target_region = 1, p1 = 0.5, p0 = 0.5)
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

  if (!is.data.frame(dat)) stop("dat must be a data frame")
  for (value in list(compute_se, region_wise)) {
    if (!is.logical(value) || length(value) != 1 || is.na(value))
      stop("compute_se and region_wise must be TRUE or FALSE")
  }
  if (!is.list(config) || (length(config) && (is.null(names(config)) ||
      any(!names(config) %in% "check_p_threshold") || anyDuplicated(names(config)))))
    stop("config may only contain check_p_threshold")

  default_config <- list(
    check_p_threshold = 0.05
  )
  config <- modifyList(default_config, config)

  required_cols <- c("y", "trt", "region")
  missing_cols  <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("dat is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (length(target_region) != 1L || is.na(target_region))
    stop("target_region must be a single nonmissing label")
  if (!(target_region %in% unique(dat$region))) {
    stop("target_region = ", target_region, " not found in dat$region")
  }

  if (length(unique(dat$region)) < 2) {
    stop("dat must contain at least 2 regions")
  }

  if (anyNA(dat) || !is.numeric(dat$y) || any(!is.finite(dat$y)) ||
      !all(dat$trt %in% c(0,1))) stop("dat must contain complete, finite outcomes and binary trt")
  for (nm in c("p1", "p0", "pi")) {
    value <- get(nm)
    if (!is.null(value)) {
      .check_number(value, nm, 0, 1, TRUE, TRUE, FALSE)
      if (!length(value) %in% c(1L, nrow(dat)))
        stop(nm, " must be scalar or have nrow(dat) entries")
    }
  }
  for (nm in c("v1", "v0")) {
    if (is.numeric(get(nm))) .check_number(get(nm), nm, 0, lower_open=TRUE)
  }
  .check_number(config$check_p_threshold, "check_p_threshold", 0, 1, TRUE, TRUE)
  for (reg in unique(dat$region)) {
    if (length(unique(dat$trt[dat$region == reg])) < 2)
      stop("Each region must contain both treatment arms")
  }
  if (anyDuplicated(names(dat))) stop("dat column names must be unique")
  all_cols   <- names(dat)
  covar_cols <- setdiff(all_cols, c("y", "trt", "region"))
  covar_cols <- covar_cols[sapply(dat[covar_cols], is.numeric)]

  if (length(covar_cols) == 0) {
    stop("No numeric covariates found in dat")
  }

  if (any(!vapply(dat[covar_cols], function(x) all(is.finite(x)), logical(1))))
    stop("Covariates must be finite")
  if (any(c("D", "D_src", "resid_sq") %in% covar_cols))
    stop("D, D_src and resid_sq are reserved internal column names")
  if (any(make.names(covar_cols) != covar_cols))
    stop("Covariate names must be syntactically valid R names")
  all_regions       <- sort(unique(dat$region))
  source_candidates <- setdiff(all_regions, target_region)

  dat$D <- as.integer(dat$region == target_region)
  q     <- mean(dat$D)

  # Fit using all candidate source regions.
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

  # Screen source regions separately in each arm.
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

  idx_wf <- which(dat$region %in% c(target_region, source_pool))
  dat_wf <- dat[idx_wf, , drop=FALSE]
  q_wf <- mean(dat_wf$D)
  subset_score <- function(x) if (is.null(x) || length(x) == 1) x else x[idx_wf]
  pi_wf <- pi
  p0_wf <- subset_score(p0)
  if (length(source_pool) == 0) {
    # Target-only AIPW: source propensity is irrelevant and region probability is one.
    pi_wf <- 1
    p0_wf <- 0.5
  } else if (length(idx_wf) < nrow(dat)) {
    # Selecting sources changes P(D=1|X); the original pi is no longer applicable.
    pi_wf <- NULL
    if (!is.null(pi)) warning("Re-estimating pi after source-region selection", call.=FALSE)
    if (is.null(p0)) p0_wf <- NULL
  }
  wf_est <- compute_dr_estimate(dat_wf, q_wf, covar_cols, compute_se=compute_se,
    p1=subset_score(p1), p0=p0_wf, pi=pi_wf, v1=v1, v0=v0)

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

#' Print a regional ATE estimate
#' @noRd
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