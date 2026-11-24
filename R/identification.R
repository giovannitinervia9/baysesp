#' @include sensitivity.R
NULL

# Partial identification with never-verified profiles ---------------------------
#
# Specification, Proposition 7.1. With P0 the set of never-verified profiles and
# R a rule that declares no profile of P0 positive:
#   TPR_R = A_R / prevalence,   FPR_R = B_R / (1 - prevalence),
#   A_R = sum over R-positive profiles of omega * psi, B_R = ... of omega * (1 - psi),
# and the prevalence is c + sum over P0 of omega * q, where c = sum over the
# other profiles of omega * psi; so it lies in [c, c + omega_0], omega_0 being
# the frequency of the profiles in P0. Under MAR psi = phi1 outside P0, and
# A_R, B_R and c are identified. Under MNAR psi also depends on beta
# (Proposition 7.2): for a fixed odds ratio the same results hold, with psi
# computed from beta.

#' Bounds for the accuracy when some profiles are never verified
#'
#' When some response profiles are never verified, the data do not identify
#' the prevalence, the sensitivities and the specificities: they depend on the
#' proportion of diseased in those profiles, which can be anything between 0
#' and 1. These are the bounds obtained by letting it vary over the whole range
#' (Proposition 7.1 of the specification). They do not depend on the prior
#' given to it.
#'
#' With [mnar()] verification the bounds also depend on the odds ratio of
#' verification between diseased and healthy subjects, which the data cannot
#' estimate either: they are computed for a fixed value, `odds_ratio`. Repeat
#' the computation over a range of values to see how they move.
#'
#' @details For a rule \eqn{R} that declares no never-verified profile positive,
#'   \deqn{A_R / (c + \omega_0) \le TPR_R \le A_R / c, \qquad
#'   B_R / (1 - c) \le FPR_R \le B_R / (1 - c - \omega_0),}
#'   and \eqn{c \le \pi \le c + \omega_0}, where \eqn{\omega_0} is the frequency of
#'   the never-verified profiles. The upper bound of the TPR is the value of the
#'   naive analysis that considers healthy every subject with a never-verified
#'   profile. The bounds are functions of identified parameters (and of the
#'   odds ratio), so they have a posterior distribution: the result gives their
#'   posterior medians and a credible region that contains the whole
#'   identification interval with probability at least `level`. For rules that
#'   declare a never-verified profile positive the bounds are `NA`.
#'
#' @inheritParams test_accuracy
#' @param fit A model fitted with [baysesp()] and [saturated()].
#' @param rules Rules or test names, see [rule_accuracy()]. By default, the
#'   individual tests.
#' @param odds_ratio Odds ratio of verification between diseased and healthy
#'   subjects, \eqn{\exp(\beta)}, see [mnar()]. By default 1 (MAR) for a model
#'   with [mar()] verification; required with [mnar()].
#' @return A data frame with columns `quantity`, `lower_bound`, `upper_bound`
#'   (posterior medians of the bounds), `lower` and `upper` (the lower quantile
#'   of the lower bound and the upper quantile of the upper bound).
#' @seealso [identified_quantities()], [sensitivity_analysis()]
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
#' identification_bounds(fit, odds_ratio = 2)
#' }
#' @export
identification_bounds <- function(fit, rules = fit@data@tests, odds_ratio = NULL,
                                  level = 0.95) {
  x <- identification_draws(fit, rules, odds_ratio, level)
  ok <- x$avoids_never
  lower <- by_rule(list(TPR = x$A / (x$c + x$omega0), FPR = x$B / (1 - x$c)), list(ok, ok))
  upper <- by_rule(list(TPR = x$A / x$c, FPR = x$B / (1 - x$c - x$omega0)), list(ok, ok))
  lower_draws <- cbind(prevalence = x$c, lower$draws)
  upper_draws <- cbind(prevalence = x$c + x$omega0, upper$draws)
  usable <- c(TRUE, lower$identified)

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
#' results do not depend on the prior on the never-verified profiles. With
#' [mnar()] verification they depend on the odds ratio of verification, and
#' are computed for a fixed value, `odds_ratio`.
#'
#' @inheritParams identification_bounds
#' @return A data frame with one row per ratio and pair of rules and columns
#'   `quantity`, `mean`, `median`, `lower`, `upper` and `probability`, the
#'   posterior probability that the ratio is greater than 1: for TPR, that the
#'   first rule is more sensitive; for FPR, that it is less specific. Only the
#'   rules that declare no never-verified profile positive are compared.
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
identified_quantities <- function(fit, rules = fit@data@tests, odds_ratio = NULL,
                                  level = 0.95) {
  x <- identification_draws(fit, rules, odds_ratio, level)
  keep <- which(x$avoids_never)
  if (length(keep) < 2L) {
    stop("At least two rules that declare no never-verified profile positive are needed.",
         call. = FALSE)
  }
  pairs <- utils::combn(keep, 2L)
  names <- colnames(x$A)
  draws <- do.call(cbind, lapply(seq_len(ncol(pairs)), function(i) {
    a <- pairs[1, i]
    b <- pairs[2, i]
    r <- cbind(x$A[, a] / x$A[, b], x$B[, a] / x$B[, b])
    colnames(r) <- sprintf("%s[%s] / %s[%s]", c("TPR", "FPR"), names[a],
                           c("TPR", "FPR"), names[b])
    r
  }))
  out <- summarise_draws(draws, level, identified = TRUE)
  out$identified <- NULL
  out$probability <- colMeans(draws > 1)
  out
}

# A_R, B_R, c and omega_0 of Proposition 7.1, draw by draw, for a fixed odds ratio.
identification_draws <- function(fit, rules, odds_ratio, level) {
  check_saturated_fit(fit, "This function")
  check_level(level)
  beta <- log(resolve_odds_ratio(fit, odds_ratio))
  rmat <- rule_matrix(as_rules(rules), fit@data)
  never <- fit@data@never_verified
  t <- transparent_draws(fit)
  psi <- profile_prevalence(t, beta, never)
  healthy <- t$omega * (1 - psi)
  healthy[, never] <- 0
  A <- (t$omega * psi) %*% rmat
  B <- healthy %*% rmat
  colnames(A) <- colnames(B) <- colnames(rmat)
  list(A = A, B = B, c = rowSums(t$omega * psi),
       omega0 = rowSums(t$omega[, never, drop = FALSE]),
       avoids_never = avoids_never(rmat, fit@data))
}

resolve_odds_ratio <- function(fit, odds_ratio) {
  if (is.null(odds_ratio)) {
    if (S7::S7_inherits(fit@model@verification, mar)) return(1)
    stop("With MNAR verification the result depends on the odds ratio of ",
         "verification: give it with `odds_ratio` (1 is MAR).", call. = FALSE)
  }
  if (!is_scalar_number(odds_ratio) || odds_ratio <= 0) {
    stop("`odds_ratio` must be a single positive number.", call. = FALSE)
  }
  odds_ratio
}
