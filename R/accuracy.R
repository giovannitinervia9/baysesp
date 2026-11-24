#' Accuracy of the individual tests
#'
#' Posterior summaries of the prevalence and of the sensitivity (`Se`),
#' specificity (`Sp`), positive predictive value (`PPV`) and negative
#' predictive value (`NPV`) of each test.
#'
#' @param fit A fitted model, from [baysesp()].
#' @param level Probability of the credible intervals.
#' @return A data frame with one row per quantity and columns `quantity`,
#'   `mean`, `median`, `lower` and `upper` (equal-tailed credible interval) and
#'   `identified`: `TRUE` if the quantity is estimated from the data alone,
#'   `FALSE` if it also depends on the prior of a parameter that the data cannot
#'   estimate (the proportion of diseased in never-verified profiles).
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- simulate_study(n = 1000, prevalence = 0.2,
#'                     sensitivity = c(A = 0.9, B = 0.8),
#'                     specificity = c(A = 0.85, B = 0.95),
#'                     verification_prob = c(0.2, 0.9, 0.9, 1))
#' d <- baysesp_data(d, tests = c("A", "B"), verified = "verified", disease = "disease")
#' test_accuracy(baysesp(d, seed = 1), level = 0.9)
#' }
#' @export
test_accuracy <- function(fit, level = 0.95) {
  if (!S7::S7_inherits(fit, baysesp_fit)) {
    stop("`fit` must be a model fitted with baysesp().", call. = FALSE)
  }
  check_level(level)
  d <- core_draws(fit)
  profiles <- fit@data@profiles
  p <- d$prevalence
  se <- d$theta_diseased %*% profiles
  sp <- 1 - d$theta_healthy %*% profiles
  ppv <- p * se / (p * se + (1 - p) * (1 - sp))
  npv <- (1 - p) * sp / ((1 - p) * sp + p * (1 - se))

  tests <- fit@data@tests
  # one block of columns per test: Se, Sp, PPV, NPV
  x <- cbind(p, do.call(cbind, lapply(seq_along(tests), function(k)
    cbind(se[, k], sp[, k], ppv[, k], npv[, k]))))
  colnames(x) <- c("prevalence",
                   paste0(c("Se", "Sp", "PPV", "NPV"), "[", rep(tests, each = 4), "]"))

  # With MAR, everything is identified unless some profile is never verified.
  # Even then the PPV of a test is identified if every profile in which the test
  # is positive has verified subjects: PPV = P(D = 1 | T_k = 1) does not involve
  # the never-verified profiles (specification, Proposition 7.1).
  mar_model <- S7::S7_inherits(fit@model@verification, mar)
  never <- fit@data@never_verified
  ppv_identified <- vapply(seq_along(tests), function(k) !any(never & profiles[, k] == 1L),
                           logical(1))
  identified <- c(!any(never), rbind(!any(never), !any(never), ppv_identified, !any(never)))

  out <- summarise_draws(x, level)
  out$identified <- mar_model & identified
  out
}

check_level <- function(level) {
  if (!is_scalar_number(level) || level <= 0 || level >= 1) {
    stop("`level` must be a single number in (0, 1).", call. = FALSE)
  }
}

# Posterior summaries of the columns of a draws x quantities matrix.
summarise_draws <- function(x, level) {
  tail <- (1 - level) / 2
  data.frame(
    quantity = colnames(x),
    mean = colMeans(x),
    median = apply(x, 2, stats::median),
    lower = apply(x, 2, stats::quantile, probs = tail, names = FALSE),
    upper = apply(x, 2, stats::quantile, probs = 1 - tail, names = FALSE),
    row.names = NULL
  )
}
