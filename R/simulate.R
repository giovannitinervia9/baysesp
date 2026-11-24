#' Simulate a study with partial verification
#'
#' Generates individual data from the selection model of the package: disease
#' status, results of `K` conditionally independent binary tests and
#' verification, which depends on the test results and, if `odds_ratio` is not
#' 1, also on the disease status. Useful for examples, tests and simulation
#' studies.
#'
#' @param n Number of subjects.
#' @param prevalence Probability of disease.
#' @param sensitivity,specificity Sensitivity and specificity of each test, one
#'   value per test. The names of `sensitivity`, if present, become the names
#'   of the test columns; otherwise the tests are called `T1`, ..., `TK`.
#' @param verification_prob Probability of verification of the healthy
#'   subjects, either a single value or one value per response profile, in
#'   binary order (see [baysesp_data()]).
#' @param odds_ratio Odds ratio of verification between diseased and healthy
#'   subjects with the same profile. The default, 1, gives verification that is
#'   missing at random; other values give verification that depends on the
#'   disease (missing not at random).
#'
#' @return A data frame with one row per subject: a 0/1 column for each test,
#'   `verified` (0/1) and `disease` (0/1, `NA` for unverified subjects).
#' @examples
#' set.seed(1)
#' d <- simulate_study(n = 500, prevalence = 0.2,
#'                     sensitivity = c(0.9, 0.7), specificity = c(0.85, 0.95),
#'                     verification_prob = c(0.1, 0.9, 0.9, 1))
#' head(d)
#' @export
simulate_study <- function(n, prevalence, sensitivity, specificity,
                           verification_prob, odds_ratio = 1) {
  is_prob <- function(x) is.numeric(x) && !anyNA(x) && all(x >= 0 & x <= 1)
  if (!is_scalar_number(n) || n < 1 || n != round(n)) {
    stop("`n` must be a positive whole number.", call. = FALSE)
  }
  if (!is_prob(prevalence) || length(prevalence) != 1L) {
    stop("`prevalence` must be a single probability.", call. = FALSE)
  }
  if (!is_prob(sensitivity) || !is_prob(specificity) || length(sensitivity) == 0L ||
      length(sensitivity) != length(specificity)) {
    stop("`sensitivity` and `specificity` must be probabilities of the same length.",
         call. = FALSE)
  }
  K <- length(sensitivity)
  J <- 2L^K
  if (!is_prob(verification_prob) || !length(verification_prob) %in% c(1L, J)) {
    stop("`verification_prob` must be a single probability or one per profile (",
         J, ").", call. = FALSE)
  }
  if (!is_scalar_number(odds_ratio) || odds_ratio <= 0) {
    stop("`odds_ratio` must be a single positive number.", call. = FALSE)
  }
  tests <- names(sensitivity)
  if (is.null(tests)) tests <- paste0("T", seq_len(K))
  if (any(c("verified", "disease") %in% tests)) {
    stop("Tests cannot be called `verified` or `disease`.", call. = FALSE)
  }

  disease <- stats::rbinom(n, 1L, prevalence)
  # probability of a positive result: sensitivity if diseased, 1 - specificity if not
  p_positive <- outer(disease, sensitivity) + outer(1L - disease, 1 - specificity)
  results <- matrix(stats::rbinom(n * K, 1L, p_positive), n, K)
  profile <- 1L + drop(results %*% 2L^(seq_len(K) - 1L))
  lambda <- rep_len(verification_prob, J)[profile]
  lambda[disease == 1L] <- stats::plogis(stats::qlogis(lambda[disease == 1L]) +
                                           log(odds_ratio))
  verified <- stats::rbinom(n, 1L, lambda)

  out <- as.data.frame(results)
  names(out) <- tests
  out$verified <- verified
  out$disease <- ifelse(verified == 1L, disease, NA_integer_)
  out
}
