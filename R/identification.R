#' @include accuracy.R
NULL

# Partial identification with never-verified profiles ---------------------------
#
# Specification, Proposition 7.1. Under MAR, with P0 the set of never-verified
# profiles and R a rule that declares no profile of P0 positive:
#   TPR_R = A_R / prevalence,   FPR_R = B_R / (1 - prevalence),
# with A_R and B_R identified. The prevalence is c + sum over P0 of omega * q,
# where c = sum over the other profiles of omega * phi1, so it lies in
# [c, c + omega_0], omega_0 being the frequency of the profiles in P0.

#' Bounds for the accuracy when some profiles are never verified
#'
#' When some response profiles are never verified, the data do not identify
#' the prevalence, the sensitivities and the specificities: they depend on the
#' proportion of diseased in those profiles, which can be anything between 0
#' and 1. These are the bounds obtained by letting it vary over the whole range
#' (Proposition 7.1 of the specification). They do not depend on the prior
#' given to it.
#'
#' @details For a rule \eqn{R} that declares no never-verified profile positive,
#'   \deqn{A_R / (c + \omega_0) \le TPR_R \le A_R / c, \qquad
#'   B_R / (1 - c) \le FPR_R \le B_R / (1 - c - \omega_0),}
#'   and \eqn{c \le \pi \le c + \omega_0}, where \eqn{\omega_0} is the frequency of
#'   the never-verified profiles. The upper bound of the TPR is the value of the
#'   naive analysis that considers healthy every subject with a never-verified
#'   profile. The bounds are functions of identified parameters, so they have a
#'   posterior distribution: the result gives their posterior medians and a
#'   credible region that contains the whole identification interval with
#'   probability at least `level`. For rules that declare a never-verified
#'   profile positive the bounds are `NA`.
#'
#' @inheritParams test_accuracy
#' @param rules Rules or test names, see [rule_accuracy()]. By default, the
#'   individual tests.
#' @return A data frame with columns `quantity`, `lower_bound`, `upper_bound`
#'   (posterior medians of the bounds), `lower` and `upper` (the lower quantile
#'   of the lower bound and the upper quantile of the upper bound).
#' @seealso [identified_quantities()]
#' @examples
#' \donttest{
#' set.seed(1)
#' # profile (0,0) is never verified
#' d <- simulate_study(n = 2000, prevalence = 0.2,
#'                     sensitivity = c(A = 0.9, B = 0.8),
#'                     specificity = c(A = 0.85, B = 0.95),
#'                     verification_prob = c(0, 0.9, 0.9, 1))
#' d <- baysesp_data(d, tests = c("A", "B"), verified = "verified", disease = "disease")
#' fit <- baysesp(d, saturated(unverified_prevalence_prior = elicit_beta(0.03, 0.06)),
#'                seed = 1)
#' identification_bounds(fit)
#' }
#' @export
identification_bounds <- function(fit, rules = fit@data@tests, level = 0.95) {
  x <- accuracy_draws(fit, rules, level)
  check_mar(fit, "identification_bounds()")
  d <- core_draws(fit)
  never <- fit@data@never_verified
  p <- x$prevalence
  c0 <- p * rowSums(d$theta_diseased[, !never, drop = FALSE])
  omega0 <- rowSums((p * d$theta_diseased + (1 - p) * d$theta_healthy)[, never, drop = FALSE])
  A <- p * x$tpr
  B <- (1 - p) * x$fpr
  ok <- avoids_never(rule_matrix(as_rules(rules), fit@data), fit@data)
  bounds <- by_rule(list(TPR = A / (c0 + omega0), FPR = B / (1 - c0)), list(ok, ok))
  upper <- by_rule(list(TPR = A / c0, FPR = B / (1 - c0 - omega0)), list(ok, ok))
  lower_draws <- cbind(prevalence = c0, bounds$draws)
  upper_draws <- cbind(prevalence = c0 + omega0, upper$draws)
  usable <- c(TRUE, bounds$identified)

  tail <- (1 - level) / 2
  na_if <- function(v) ifelse(usable, v, NA_real_)
  data.frame(
    quantity = colnames(lower_draws),
    lower_bound = na_if(apply(lower_draws, 2, stats::median)),
    upper_bound = na_if(apply(upper_draws, 2, stats::median)),
    lower = na_if(apply(lower_draws, 2, stats::quantile, probs = tail, names = FALSE)),
    upper = na_if(apply(upper_draws, 2, stats::quantile, probs = 1 - tail, names = FALSE)),
    row.names = NULL
  )
}

#' Quantities identified when some profiles are never verified
#'
#' When some response profiles are never verified, the absolute accuracy of
#' the tests is not identified, but the ratios of the true positive rates and
#' of the false positive rates of two rules are, if neither rule declares a
#' never-verified profile positive (Proposition 7.1 of the specification). So
#' are the signs of the differences in sensitivity and in specificity. These
#' results do not depend on the prior on the never-verified profiles.
#'
#' @inheritParams identification_bounds
#' @return A data frame with one row per ratio and pair of rules, columns as in
#'   [compare_tests()]: `probability` is the posterior probability that the
#'   ratio is greater than 1, that is, for TPR, that the first rule is more
#'   sensitive, and for FPR, that it is less specific. Only the rules that
#'   declare no never-verified profile positive are compared.
#' @seealso [identification_bounds()], [compare_tests()]
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- simulate_study(n = 2000, prevalence = 0.2,
#'                     sensitivity = c(A = 0.9, B = 0.8),
#'                     specificity = c(A = 0.85, B = 0.95),
#'                     verification_prob = c(0, 0.9, 0.9, 1))
#' d <- baysesp_data(d, tests = c("A", "B"), verified = "verified", disease = "disease")
#' fit <- baysesp(d, saturated(unverified_prevalence_prior = elicit_beta(0.03, 0.06)),
#'                seed = 1)
#' identified_quantities(fit, rules = list("A", "B", any_positive()))
#' }
#' @export
identified_quantities <- function(fit, rules = fit@data@tests, level = 0.95) {
  x <- accuracy_draws(fit, rules, level)
  check_mar(fit, "identified_quantities()")
  keep <- which(x$ppv_identified)
  if (length(keep) < 2L) {
    stop("At least two rules that declare no never-verified profile positive are needed.",
         call. = FALSE)
  }
  pairs <- utils::combn(keep, 2L)
  out <- do.call(rbind, lapply(seq_len(ncol(pairs)), function(i) {
    r <- pairs[, i]
    rules <- list(as_rules(rules)[[r[1]]], as_rules(rules)[[r[2]]])
    compare_tests(fit, rules[[1]], rules[[2]], level)[4:5, ]
  }))
  out$identified <- NULL
  rownames(out) <- NULL
  out
}

check_mar <- function(fit, what) {
  if (!S7::S7_inherits(fit@model@verification, mar)) {
    stop(what, " assumes MAR verification.", call. = FALSE)
  }
}
