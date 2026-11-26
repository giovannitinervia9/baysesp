#' @include rules.R fit.R
NULL

# Accuracy of tests and rules ---------------------------------------------------
#
# Every quantity is a function of the prevalence and of the distributions of the
# profiles in the diseased and in the healthy (specification, section 4): for a
# rule with vector r, TPR = r' theta_1 and FPR = r' theta_0. Each function
# computes them draw by draw and summarises the posterior.
#
# Identification (specification, Proposition 7.1), with MAR verification:
# - with no never-verified profile, everything is identified;
# - otherwise the absolute quantities are not, but for rules that declare no
#   never-verified profile positive the PPV and the ratios of TPR and FPR
#   between two such rules are.

#' Accuracy of the individual tests
#'
#' Posterior summaries of the prevalence and of the sensitivity (`Se`),
#' specificity (`Sp`), positive predictive value (`PPV`) and negative
#' predictive value (`NPV`) of each test.
#'
#' @template accuracy
#' @details The first row is the prevalence \eqn{\pi}, the proportion of
#'   diseased subjects in the population.
#' @param fit A fitted model, from [baysesp()].
#' @param level Probability of the credible intervals.
#' @return A data frame with one row per quantity and columns `quantity`,
#'   `mean`, `median`, `lower` and `upper` (equal-tailed credible interval) and
#'   `identified`: `TRUE` if the quantity is estimated from the data alone,
#'   `FALSE` if it also depends on the prior of a parameter that the data cannot
#'   estimate (the proportion of diseased in never-verified profiles).
#' @seealso [rule_accuracy()] for combinations of tests, [compare_tests()].
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
  x <- accuracy_draws(fit, fit@data@tests, level)
  out <- by_rule(list(Se = x$tpr, Sp = 1 - x$fpr, PPV = x$ppv, NPV = x$npv),
                 list(x$all, x$all, x$ppv_identified, x$all))
  summarise_draws(cbind(prevalence = x$prevalence, out$draws), level,
                  identified = c(x$all[1], out$identified))
}

#' Accuracy of diagnostic rules
#'
#' Posterior summaries of the sensitivity (`Se`), specificity (`Sp`),
#' positive predictive value (`PPV`) and negative predictive value (`NPV`) of
#' rules that combine the tests. For a single test the result is the same as
#' in [test_accuracy()].
#'
#' @template accuracy
#' @inheritParams test_accuracy
#' @param rules A rule, such as [any_positive()], the name of a test, or a list
#'   of them.
#' @return A data frame as in [test_accuracy()], without the row of the
#'   prevalence.
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- simulate_study(n = 1000, prevalence = 0.2,
#'                     sensitivity = c(A = 0.9, B = 0.8),
#'                     specificity = c(A = 0.85, B = 0.95),
#'                     verification_prob = c(0.2, 0.9, 0.9, 1))
#' d <- baysesp_data(d, tests = c("A", "B"), verified = "verified", disease = "disease")
#' fit <- baysesp(d, seed = 1)
#' rule_accuracy(fit, list(any_positive(), all_positive()))
#' }
#' @export
rule_accuracy <- function(fit, rules, level = 0.95) {
  x <- accuracy_draws(fit, rules, level)
  out <- by_rule(list(Se = x$tpr, Sp = 1 - x$fpr, PPV = x$ppv, NPV = x$npv),
                 list(x$all, x$all, x$ppv_identified, x$all))
  summarise_draws(out$draws, level, identified = out$identified)
}

#' Compare two tests or rules
#'
#' Posterior summaries of the differences in sensitivity, specificity and
#' Youden index (`first` minus `second`) and of the ratios of the true and false
#' positive rates (`first` over `second`).
#'
#' @template accuracy
#' @details For two rules \eqn{R} (`first`) and \eqn{R'} (`second`):
#'
#' - **differences**: \eqn{\Delta Se = Se_R - Se_{R'}},
#'   \eqn{\Delta Sp = Sp_R - Sp_{R'}} and, for the Youden index
#'   \eqn{J_R = Se_R + Sp_R - 1}, which summarises the accuracy in one number
#'   (0 for a useless test, 1 for a perfect one),
#'   \eqn{\Delta J = \Delta Se + \Delta Sp}. A positive difference favours the
#'   first rule;
#' - **ratios** of the true positive rates,
#'   \eqn{rTPR = Se_R / Se_{R'}}, and of the false positive rates,
#'   \eqn{rFPR = (1 - Sp_R) / (1 - Sp_{R'})}. For example, \eqn{rTPR = 1.2} means
#'   that the first rule detects 20% more diseased subjects; \eqn{rFPR < 1}
#'   means that it gives fewer false positives.
#'
#' The column `probability` is the posterior probability that a difference is
#' positive, or that a ratio is greater than 1. \eqn{P(rTPR > 1)} equals
#' \eqn{P(\Delta Se > 0)}, and \eqn{P(rFPR > 1)} equals \eqn{P(\Delta Sp < 0)}.
#'
#' When some profile is never verified, the differences depend on the prior of
#' its proportion of diseased, but the ratios do not, as long as neither rule
#' declares a never-verified profile positive. In that case the signs of the differences in sensitivity
#' and in specificity are also identified, through the probabilities of the
#' ratios.
#' @inheritParams test_accuracy
#' @param first,second Two tests (by name) or two rules, see [test_rule].
#' @return A data frame as in [test_accuracy()], with the additional column
#'   `probability`: the posterior probability that a difference is positive or
#'   that a ratio is greater than 1.
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- simulate_study(n = 1000, prevalence = 0.2,
#'                     sensitivity = c(A = 0.9, B = 0.8),
#'                     specificity = c(A = 0.85, B = 0.95),
#'                     verification_prob = c(0.2, 0.9, 0.9, 1))
#' d <- baysesp_data(d, tests = c("A", "B"), verified = "verified", disease = "disease")
#' fit <- baysesp(d, seed = 1)
#' compare_tests(fit, "A", "B")
#' compare_tests(fit, any_positive(), "A")
#' }
#' @export
compare_tests <- function(fit, first, second, level = 0.95) {
  rules <- c(as_rules(first, "first"), as_rules(second, "second"))
  if (length(rules) != 2L) {
    stop("`first` and `second` must be one test or rule each.", call. = FALSE)
  }
  x <- accuracy_draws(fit, rules, level)
  tpr <- x$tpr
  fpr <- x$fpr
  diff_se <- tpr[, 1] - tpr[, 2]
  diff_sp <- fpr[, 2] - fpr[, 1]
  draws <- cbind(diff_se, diff_sp, diff_se + diff_sp, tpr[, 1] / tpr[, 2], fpr[, 1] / fpr[, 2])
  a <- colnames(tpr)[1]
  b <- colnames(tpr)[2]
  colnames(draws) <- c(sprintf("%s[%s] - %s[%s]", c("Se", "Sp", "Youden"), a,
                               c("Se", "Sp", "Youden"), b),
                       sprintf("%s[%s] / %s[%s]", c("TPR", "FPR"), a, c("TPR", "FPR"), b))
  ratios_identified <- all(x$ppv_identified)
  out <- summarise_draws(draws, level,
                         identified = c(rep(x$all[1], 3), rep(ratios_identified, 2)))
  out$probability <- colMeans(sweep(draws, 2, c(0, 0, 0, 1, 1), `>`))
  out
}

# Draws of prevalence, TPR, FPR, PPV and NPV of a set of rules, with the
# identification of each rule.
accuracy_draws <- function(fit, rules, level) {
  if (!S7::S7_inherits(fit, baysesp_fit)) {
    stop("`fit` must be a model fitted with baysesp().", call. = FALSE)
  }
  check_level(level)
  rmat <- rule_matrix(as_rules(rules), fit@data)
  d <- core_draws(fit)
  p <- d$prevalence
  tpr <- d$theta_diseased %*% rmat
  fpr <- d$theta_healthy %*% rmat
  colnames(tpr) <- colnames(fpr) <- colnames(rmat)
  mar_model <- S7::S7_inherits(fit@model@verification, mar)
  list(
    prevalence = p, tpr = tpr, fpr = fpr,
    ppv = p * tpr / (p * tpr + (1 - p) * fpr),
    npv = (1 - p) * (1 - fpr) / ((1 - p) * (1 - fpr) + p * (1 - tpr)),
    all = rep(mar_model && !any(fit@data@never_verified), ncol(rmat)),
    ppv_identified = mar_model & avoids_never(rmat, fit@data)
  )
}

# Rules that declare no never-verified profile positive.
avoids_never <- function(rmat, data) {
  colSums(rmat[data@never_verified, , drop = FALSE]) == 0
}

# Binds draws x rules matrices, one per quantity, into one draws matrix with the
# columns grouped by rule ("Se[A]", "Sp[A]", ..., "Se[B]", ...).
by_rule <- function(blocks, identified) {
  B <- length(blocks)
  R <- ncol(blocks[[1]])
  order <- as.vector(outer((seq_len(B) - 1L) * R, seq_len(R), `+`))
  draws <- do.call(cbind, unname(blocks))[, order, drop = FALSE]
  colnames(draws) <- paste0(names(blocks), "[", rep(colnames(blocks[[1]]), each = B), "]")
  list(draws = draws, identified = unlist(identified)[order])
}

check_level <- function(level) {
  if (!is_scalar_number(level) || level <= 0 || level >= 1) {
    stop("`level` must be a single number in (0, 1).", call. = FALSE)
  }
}

# Posterior summaries of the columns of a draws x quantities matrix.
summarise_draws <- function(x, level, identified) {
  tail <- (1 - level) / 2
  data.frame(
    quantity = colnames(x),
    mean = colMeans(x),
    median = apply(x, 2, stats::median),
    lower = apply(x, 2, stats::quantile, probs = tail, names = FALSE),
    upper = apply(x, 2, stats::quantile, probs = 1 - tail, names = FALSE),
    identified = identified,
    row.names = NULL
  )
}
