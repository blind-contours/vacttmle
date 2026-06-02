###############################################################################
# Shared utilities for the target-aligned VA-CT-TMLE simulation study.
###############################################################################

`%||%` <- function(x, y) if (is.null(x)) y else x

va_expit <- function(x) plogis(pmax(pmin(x, 35), -35))

va_logit <- function(p, eps = 1e-8) {
  p <- pmin(pmax(p, eps), 1 - eps)
  qlogis(p)
}

va_bound <- function(x, eps = 1e-8) pmin(pmax(x, eps), 1 - eps)

va_project_dir <- function() {
  getOption("va_ct_tmle.project_dir", default = getwd())
}

va_path <- function(...) file.path(va_project_dir(), ...)

va_ensure_output_dirs <- function() {
  dirs <- c(
    "outputs/toy_checks", "outputs/truth", "outputs/pilot",
    "outputs/main", "outputs/gammaZ", "outputs/jointZ",
    "outputs/jointZ_adaptive", "outputs/comparators", "outputs/s6",
    "outputs/tables", "outputs/figures"
  )
  invisible(lapply(va_path(dirs), dir.create, recursive = TRUE, showWarnings = FALSE))
}

va_package_versions <- function() {
  pkgs <- c("data.table", "ggplot2", "testthat", "knitr", "survival")
  data.table(
    package = pkgs,
    version = vapply(pkgs, function(p) {
      if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_
    }, character(1))
  )
}

va_build_formula <- function(outcome, data, candidates, offset_var = NULL) {
  keep <- character()
  for (v in candidates) {
    if (v == "k_f") {
      if ("k_f" %in% names(data) && uniqueN(data$k_f) > 1) keep <- c(keep, "k_f")
    } else if (v %in% names(data)) {
      z <- data[[v]]
      if (sum(!is.na(z)) > 1 && uniqueN(z[!is.na(z)]) > 1) keep <- c(keep, v)
    }
  }
  rhs <- if (length(keep)) paste(keep, collapse = " + ") else "1"
  if (!is.null(offset_var)) rhs <- paste(rhs, "+ offset(", offset_var, ")", sep = "")
  as.formula(paste(outcome, "~", rhs))
}

va_fit_poisson_rate <- function(data, outcome, candidates) {
  d <- copy(data)
  d <- d[time_at_risk > 1e-10]
  d[, log_pt := log(time_at_risk)]
  d[, k_f := factor(k)]
  total_time <- sum(d$time_at_risk)
  fallback <- (sum(d[[outcome]], na.rm = TRUE) + 0.5) / (total_time + 1)
  if (!nrow(d)) {
    return(list(model = NULL, fallback = fallback, formula = NULL))
  }
  f <- va_build_formula(outcome, d, candidates, offset_var = "log_pt")
  fit <- tryCatch(
    suppressWarnings(glm(f, family = poisson(), data = d)),
    error = function(e) NULL
  )
  list(model = fit, fallback = fallback, formula = f)
}

va_predict_rate <- function(fit, newdata) {
  nd <- copy(newdata)
  nd[, log_pt := 0]
  nd[, k_f := factor(k)]
  if (is.null(fit$model)) return(rep(fit$fallback, nrow(nd)))
  pred <- tryCatch(
    as.numeric(predict(fit$model, newdata = nd, type = "response")),
    error = function(e) rep(fit$fallback, nrow(nd))
  )
  pmax(pred, 1e-10)
}

va_solve_poisson_eps <- function(h, count, mu0) {
  ok <- is.finite(h) & is.finite(count) & is.finite(mu0) & mu0 > 0
  h <- h[ok]
  count <- count[ok]
  mu0 <- mu0[ok]
  if (!length(h) || max(abs(h)) < 1e-14) return(0)
  score <- function(eps) sum(h * (count - mu0 * exp(eps * h)))
  s0 <- score(0)
  if (!is.finite(s0) || abs(s0) < 1e-10) return(0)
  lo <- -25
  hi <- 25
  slo <- score(lo)
  shi <- score(hi)
  if (is.finite(slo) && is.finite(shi) && slo * shi <= 0) {
    return(uniroot(score, c(lo, hi), tol = 1e-8)$root)
  }
  opt <- optimize(function(eps) -sum(count * eps * h - mu0 * exp(eps * h)),
                  interval = c(lo, hi))
  opt$minimum
}

va_solve_logistic_eps <- function(offset, h, y) {
  ok <- is.finite(offset) & is.finite(h) & is.finite(y)
  offset <- offset[ok]
  h <- h[ok]
  y <- y[ok]
  y <- va_bound(y)
  if (!length(h) || max(abs(h)) < 1e-14) return(0)
  score <- function(eps) sum(h * (y - va_expit(offset + eps * h)))
  s0 <- score(0)
  if (!is.finite(s0) || abs(s0) < 1e-10) return(0)
  lo <- -50
  hi <- 50
  slo <- score(lo)
  shi <- score(hi)
  if (is.finite(slo) && is.finite(shi) && slo * shi <= 0) {
    return(uniroot(score, c(lo, hi), tol = 1e-8)$root)
  }
  opt <- optimize(function(eps) {
    p <- va_bound(va_expit(offset + eps * h))
    -sum(y * log(p) + (1 - y) * log1p(-p))
  }, interval = c(lo, hi))
  opt$minimum
}

va_effective_sample_size <- function(w) {
  if (!length(w) || sum(w, na.rm = TRUE) <= 0) return(NA_real_)
  sum(w, na.rm = TRUE)^2 / sum(w^2, na.rm = TRUE)
}

va_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(as.data.table(x), path)
  invisible(path)
}
