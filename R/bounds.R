#' @include sensitivity.R
NULL

# Bounds and ratios with never-verified profiles ---------------------------
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

#' Bounds for the absolute accuracy when some profiles are never verified
#'
#' Gives the range of values that the prevalence and the sensitivity and
#' specificity of each test (or rule) can take, whatever the proportion of
#' diseased subjects in the response profiles that are never verified.
#'
#' When some profiles (combinations of test results) are never verified, the
#' data say nothing about how many diseased subjects they contain, and the
#' absolute accuracy depends on it. A prior on that proportion gives a single
#' answer, which rests entirely on the prior (Buzoianu and Kadane, 2008). These bounds
#' instead show what the data alone allow (Gustafson, 2015): the lowest and highest value of each quantity when the
#' proportion of diseased in the never-verified profiles ranges from 0 to 1.
#' They do not depend on the prior. Wide bounds mean that the absolute accuracy
#' is mostly determined by the assumption on the never-verified profiles.
#'
#' The bounds also depend on the odds ratio of verification between diseased
#' and healthy subjects with the same profile, which the data cannot estimate
#' either: they are computed for one fixed value, `odds_ratio` (by default 1,
#' verification that does not depend on the disease). Call the function again
#' with other values to see how the bounds move.
#'
#' @template accuracy
#' @details **Bounds**. Let
#'   \eqn{\psi_t = P(D = 1 \mid T = t)} be the proportion of diseased in profile
#'   \eqn{t} and \eqn{\omega_t = P(T = t)} its frequency. For a profile with
#'   verified subjects \eqn{\psi_t} is estimated from the data: with odds
#'   ratio of verification \eqn{OR = \exp(\beta)},
#'   \eqn{\psi_t = v_t \phi_t + (1 - v_t)\,\mathrm{expit}(\mathrm{logit}\,\phi_t - \beta)},
#'   where \eqn{v_t} is the fraction of verified subjects in the profile and
#'   \eqn{\phi_t} the proportion of diseased among them (with \eqn{OR = 1},
#'   \eqn{\psi_t = \phi_t}). For a never-verified profile \eqn{\psi_t} is
#'   unknown. For a rule \eqn{R} that declares no never-verified
#'   profile positive, let
#'   \deqn{A_R = \sum_{t: r_R(t) = 1} \omega_t \psi_t, \qquad
#'   B_R = \sum_{t: r_R(t) = 1} \omega_t (1 - \psi_t),}
#'   the proportions of subjects who are declared positive and are diseased or
#'   healthy, \eqn{c} the proportion of diseased subjects with a verified
#'   profile and \eqn{\omega_0} the frequency of the never-verified profiles.
#'   Then \eqn{Se_R = A_R / \pi} and \eqn{1 - Sp_R = B_R / (1 - \pi)}, and since
#'   \eqn{c \le \pi \le c + \omega_0},
#'   \deqn{\frac{A_R}{c + \omega_0} \le Se_R \le \frac{A_R}{c}, \qquad
#'   1 - \frac{B_R}{1 - c - \omega_0} \le Sp_R \le 1 - \frac{B_R}{1 - c}.}
#'   The upper bound of the sensitivity and the lower bound of the specificity
#'   are the values of the naive analysis that considers healthy every subject
#'   with a never-verified profile; the other bounds correspond to considering
#'   all of them diseased.
#'
#'   The bounds are functions of identified parameters (and of the odds ratio),
#'   so they have a posterior distribution: the result gives their posterior
#'   medians and a credible region that contains the whole identification
#'   interval with probability at least `level` (Gustafson, 2015). For rules that declare a
#'   never-verified profile positive the bounds are `NA`.
#' @inheritParams test_accuracy
#' @param fit A model fitted with [baysesp()] and [saturated()].
#' @param rules Rules or test names, see [rule_accuracy()]. By default, the
#'   individual tests.
#' @param odds_ratio The odds ratio of verification between diseased and
#'   healthy subjects with the same response profile, \eqn{\exp(\beta)}: a
#'   single positive number (one value only; to compare several values, call
#'   the function once for each). 1 means that verification does not depend on
#'   the disease status (missing at random); values above 1 mean that the
#'   diseased are verified more often. The default, `NULL`, uses 1 for a model
#'   fitted with `verification = mar()`, and gives an error for a model fitted
#'   with `verification = mnar()`, for which a value must be given.
#' @return A data frame with columns `quantity`, `lower_bound`, `upper_bound`
#'   (posterior medians of the bounds), `lower` and `upper` (the lower quantile
#'   of the lower bound and the upper quantile of the upper bound).
#' @references
#' Buzoianu, M. and Kadane, J. B. (2008). Adjusting for verification bias in diagnostic test evaluation: A Bayesian approach. *Statistics in Medicine*, 27(13), 2453–2473. \doi{10.1002/sim.3099}
#'
#' Gustafson, P. (2015). *Bayesian Inference for Partially Identified Models: Exploring the Limits of Limited Data*. Chapman & Hall/CRC, Boca Raton, FL.
#'
#' Torrance-Rynard, V. L. and Walter, S. D. (1997). Effects of dependent errors in the assessment of diagnostic test performance. *Statistics in Medicine*, 16(19), 2157–2175. \doi{10.1002/(SICI)1097-0258(19971015)16:19<2157::AID-SIM653>3.0.CO;2-X}
#'
#' Vacek, P. M. (1985). The effect of conditional dependence on the evaluation of diagnostic tests. *Biometrics*, 41(4), 959–968. \doi{10.2307/2530967}
#' @seealso [accuracy_ratios()], [sensitivity_analysis()]
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
#' accuracy_bounds(fit)
#' accuracy_bounds(fit, odds_ratio = 2)
#' }
#' @export
accuracy_bounds <- function(fit, rules = fit@data@tests, odds_ratio = NULL,
                                  level = 0.95) {
  x <- identification_draws(fit, rules, odds_ratio, level)
  ok <- x$avoids_never
  # Se = TPR; Sp = 1 - FPR, so its lower bound comes from the upper bound of FPR
  lower <- by_rule(list(Se = x$A / (x$c + x$omega0), Sp = 1 - x$B / (1 - x$c - x$omega0)),
                   list(ok, ok))
  upper <- by_rule(list(Se = x$A / x$c, Sp = 1 - x$B / (1 - x$c)), list(ok, ok))
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

#' Ratios between the accuracy of tests when some profiles are never verified
#'
#' Gives, for each pair of tests (or rules), the ratio of their sensitivities
#' and the ratio of their false positive rates (one minus the specificity).
#' Unlike the absolute accuracy, these ratios are estimated by the data even
#' when some response profiles are never verified (Alonzo, 2005; Pepe and Alonzo,
#' 2001; Schatzkin et al., 1987).
#'
#' When some profiles (combinations of test results) are never verified, the
#' data say nothing about how many diseased subjects they contain (Gustafson, 2015), and the
#' sensitivities and specificities depend on it. If neither test declares a
#' never-verified profile positive, however, the unknown diseased subjects are
#' missed by both tests: they lower the two sensitivities by the same factor,
#' and the ratio of the sensitivities does not change. The same holds for the
#' false positive rates. These ratios therefore tell how the tests compare,
#' for example that one detects 15% more diseased subjects than the other,
#' without any assumption on the never-verified profiles; the posterior
#' probability that a ratio is greater than 1 tells which test is more
#' sensitive, or less specific. Tests or rules that declare a never-verified
#' profile positive are left out.
#'
#' The ratios depend on the odds ratio of verification between diseased and
#' healthy subjects with the same profile, which the data cannot estimate:
#' they are computed for one fixed value, `odds_ratio` (by default 1,
#' verification that does not depend on the disease).
#'
#' @template accuracy
#' @details **Why the ratios are estimated by the data.** Let
#'   \eqn{\omega_t} be the frequency of profile \eqn{t}, \eqn{\psi_t} the
#'   proportion of diseased in it, known from the data for the profiles with
#'   verified subjects, and \eqn{\pi} the prevalence, which depends on the
#'   unknown \eqn{\psi_t} of the never-verified profiles. For a rule \eqn{R}
#'   that declares no never-verified profile positive,
#'   \deqn{A_R = \sum_{t: r_R(t) = 1} \omega_t \psi_t, \qquad
#'   B_R = \sum_{t: r_R(t) = 1} \omega_t (1 - \psi_t)}
#'   are known, and \eqn{Se_R = A_R / \pi}, \eqn{1 - Sp_R = B_R / (1 - \pi)}.
#'   For two such rules the unknown prevalence cancels out in the ratios
#'   \deqn{rTPR = \frac{Se_R}{Se_{R'}} = \frac{A_R}{A_{R'}}, \qquad
#'   rFPR = \frac{1 - Sp_R}{1 - Sp_{R'}} = \frac{B_R}{B_{R'}},}
#'   which are therefore estimated from the data alone (Alonzo, 2005; Pepe and Alonzo, 2001; Schatzkin
#'   et al., 1987). \eqn{rTPR > 1} means
#'   that \eqn{R} is more sensitive than \eqn{R'}, \eqn{rFPR > 1} that it is
#'   less specific, whatever the proportion of diseased in the never-verified
#'   profiles.
#' @inheritParams accuracy_bounds
#' @return A data frame with one row per ratio and pair of rules and columns
#'   `quantity`, `mean`, `median`, `lower`, `upper` and `probability`, the
#'   posterior probability that the ratio is greater than 1: for rTPR, that the
#'   first rule is more sensitive; for rFPR, that it is less specific. Only the
#'   rules that declare no never-verified profile positive are compared.
#' @references
#' Alonzo, T. A. (2005). Verification bias-corrected estimators of the relative true and false positive rates of two binary screening tests. *Statistics in Medicine*, 24(3), 403–417. \doi{10.1002/sim.1959}
#'
#' Gustafson, P. (2015). *Bayesian Inference for Partially Identified Models: Exploring the Limits of Limited Data*. Chapman & Hall/CRC, Boca Raton, FL.
#'
#' Pepe, M. S. and Alonzo, T. A. (2001). Comparing disease screening tests when true disease status is ascertained only for screen positives. *Biostatistics*, 2(3), 249–260. \doi{10.1093/biostatistics/2.3.249}
#'
#' Schatzkin, A., Connor, R. J., Taylor, P. R. and Bunnag, B. (1987). Comparing new and old screening tests when a reference procedure cannot be performed on all screenees. *American Journal of Epidemiology*, 125(4), 672–678. \doi{10.1093/oxfordjournals.aje.a114580}
#'
#' Torrance-Rynard, V. L. and Walter, S. D. (1997). Effects of dependent errors in the assessment of diagnostic test performance. *Statistics in Medicine*, 16(19), 2157–2175. \doi{10.1002/(SICI)1097-0258(19971015)16:19<2157::AID-SIM653>3.0.CO;2-X}
#'
#' Vacek, P. M. (1985). The effect of conditional dependence on the evaluation of diagnostic tests. *Biometrics*, 41(4), 959–968. \doi{10.2307/2530967}
#' @seealso [accuracy_bounds()], [compare_tests()]
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
#' accuracy_ratios(fit, rules = list("A", "B", any_positive()))
#' }
#' @export
accuracy_ratios <- function(fit, rules = fit@data@tests, odds_ratio = NULL,
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
