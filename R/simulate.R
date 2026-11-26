#' Simulate a study with partial verification
#'
#' Generates individual data for a study in which `K` binary tests are applied
#' to all subjects, but only some of them are verified by the reference
#' standard. For each subject the disease status is drawn first, then the `K`
#' test results and finally the verification, which depends on the response
#' profile (the vector of test results) and, when `odds_ratio` is not 1, also
#' on the disease status. Useful for examples, tests and simulation studies.
#'
#' @details
#' For each subject the function simulates:
#'
#' - the disease status \eqn{D}, Bernoulli with probability `prevalence`;
#' - the `K` test results \eqn{T_1, \dots, T_K} given \eqn{D}, with
#'   \eqn{P(T_k = 1 \mid D = 1)} equal to `sensitivity[k]` and
#'   \eqn{P(T_k = 0 \mid D = 0)} equal to `specificity[k]`. The tests are
#'   independent given \eqn{D} with `dependence = independence()`, and
#'   dependent given \eqn{D} with the other values of `dependence` (see the
#'   paragraph on the dependence between tests below);
#' - the verification \eqn{V}, Bernoulli with probability
#'   \eqn{\lambda_{t,d} = P(V = 1 \mid T = t, D = d)}.
#'
#' The verification probability of a healthy subject with profile \eqn{t} is
#' \eqn{\lambda_{t,0}} = `verification_prob` for that profile. The verification
#' probability of a diseased subject with the same profile is
#' \eqn{\lambda_{t,1} = \mathrm{expit}(\mathrm{logit}\,\lambda_{t,0} +
#' \log(\mathrm{odds\_ratio}))}, so that the odds of verification of diseased
#' subjects are `odds_ratio` times the odds for healthy subjects with the same
#' profile.
#'
#' **Dependence between the tests.** The argument `dependence` sets the joint
#' distribution of the `K` results within each disease class, that is the
#' probabilities \eqn{\theta_{d,t} = P(T = t \mid D = d)} of the response
#' profiles \eqn{t} among the diseased (\eqn{d = 1}) and among the healthy
#' (\eqn{d = 0}). Every mechanism keeps the sensitivities and specificities
#' given in `sensitivity` and `specificity` as the marginal probabilities, and
#' changes only how often the tests are right or wrong together:
#'
#' - `dependence = independence()`:
#'   \eqn{\theta_{1,t} = \prod_k Se_k^{t_k}(1 - Se_k)^{1 - t_k}} and
#'   \eqn{\theta_{0,t} = \prod_k (1 - Sp_k)^{t_k} Sp_k^{1 - t_k}};
#' - `dependence = latent_dependence(diseased_loadings, healthy_loadings)`:
#'   each subject has a latent variable \eqn{Z \sim N(0, 1)}, independent of
#'   \eqn{D}, and given \eqn{D = d} and \eqn{Z = z} the tests are independent
#'   with \eqn{P(T_k = 1 \mid D = d, Z = z) = \Phi(a_{dk} + b_{dk} z)}, where
#'   \eqn{b_{dk}} is the loading of test \eqn{k} in class \eqn{d} and
#'   \eqn{a_{dk}} is chosen so that the marginal probabilities are the
#'   sensitivities and specificities. Tests \eqn{j} and \eqn{k} are positively
#'   dependent when \eqn{b_{dj} b_{dk} > 0} and negatively dependent when
#'   \eqn{b_{dj} b_{dk} < 0};
#' - `dependence = probit_dependence(diseased_correlation, healthy_correlation)`:
#'   given \eqn{D = d}, test \eqn{k} is positive when a standard normal latent
#'   variable \eqn{Y_k} falls below a threshold that gives the required
#'   sensitivity or specificity, and the latent variables of the `K` tests are
#'   jointly normal with correlation matrix \eqn{R_d}. The element
#'   \eqn{R_{d,jk}} is the correlation, of any sign, between tests \eqn{j} and
#'   \eqn{k} in class \eqn{d}.
#'
#' The size of the loadings or of the correlation matrices is checked here
#' against the number of tests `K`, which the dependence mechanism does not
#' know when it is created. A correlation matrix that is not positive definite
#' for `K` tests gives an error.
#'
#' **Supported data-generating mechanisms.**
#'
#' - any number `K >= 1` of tests, with arbitrary prevalence, sensitivities and
#'   specificities, with tests that are independent given the disease status
#'   or dependent through one of the mechanisms above;
#' - verification that depends on the response profile only, with the same
#'   probability for all profiles (a single `verification_prob`) or with one
#'   probability per profile. A value 0 means that the profile is never
#'   verified, a value 1 that it is always verified (for healthy subjects; see
#'   `odds_ratio` for diseased subjects);
#' - missing at random verification (`odds_ratio = 1`), where verification
#'   depends only on the profile, and missing not at random verification
#'   (`odds_ratio` different from 1), where the odds of verification also
#'   depend on the disease status through the same odds ratio for all profiles.
#'   Profiles with verification probability 0 or 1 are never verified, or
#'   always verified, whatever `odds_ratio` is.
#'
#' **Not supported:** dependence between the tests given the disease status
#' other than through the two latent mechanisms above (for example log-linear
#' interactions given as odds ratios), covariates, verification that depends on something other
#' than the response profile and the disease status, and a disease-dependent
#' verification whose odds ratio differs across profiles.
#'
#' @param n Number of subjects: a single whole number, at least 1. No default.
#' @param prevalence Probability of disease: a single number in
#'   \eqn{[0, 1]}. No default.
#' @param sensitivity Sensitivity of each test: a numeric vector with values in
#'   \eqn{[0, 1]}, one element per test, of length `K >= 1`. No default. The
#'   names of `sensitivity`, if present, become the names of the test columns
#'   of the output; if it has no names, the tests are called `T1`, ..., `TK`.
#'   The names `verified` and `disease` are not allowed.
#' @param specificity Specificity of each test: a numeric vector with values in
#'   \eqn{[0, 1]} and the same length as `sensitivity`; its k-th element is the
#'   specificity of the k-th test. Its names are ignored. No default.
#' @param verification_prob Probability of verification of healthy subjects,
#'   with values in \eqn{[0, 1]}. It is either a single number, which applies to
#'   all response profiles, or a vector of length \eqn{2^K} with one number per
#'   response profile, in binary order: the profile \eqn{t} is at position
#'   \eqn{1 + \sum_k t_k 2^{k-1}}, so the first profile is the one with all tests
#'   negative and the first test changes fastest. With `K = 2` the order is
#'   (0,0), (1,0), (0,1), (1,1). A value 0 means that the profile is never
#'   verified and a value 1 that it is always verified. No default.
#' @param odds_ratio Odds ratio of verification between diseased and healthy
#'   subjects with the same response profile: a single positive number, the
#'   same for all profiles. Default `1`, which gives verification that is
#'   missing at random (it does not depend on the disease status). Values
#'   greater than 1 make diseased subjects more likely to be verified, values
#'   between 0 and 1 less likely; any value different from 1 gives verification
#'   that is missing not at random.
#' @param dependence How the test results depend on each other given the
#'   disease status: a single dependence mechanism, one of
#'   `dependence = independence()` (the tests are independent given the
#'   disease status), `dependence = latent_dependence(...)` (the tests share a
#'   normal latent variable, with one loading per test and disease class) or
#'   `dependence = probit_dependence(...)` (the tests have jointly normal
#'   latent variables, with one correlation per pair of tests and disease
#'   class). Default `independence()`. Whatever the mechanism, `sensitivity`
#'   and `specificity` remain the marginal sensitivity and specificity of each
#'   test.
#'
#' @return A data frame with one row per subject and the following columns: one
#'   0/1 column for each test, always observed, named as described in
#'   `sensitivity`; `verified` (0/1); `disease` (0/1 for verified subjects,
#'   `NA` for unverified subjects).
#' @examples
#' set.seed(1)
#'
#' # Two tests, verification depending only on the results (MAR): subjects
#' # negative to both tests are verified with probability 0.1, the others with
#' # probability 0.9 or 1. The profiles are in the order (0,0), (1,0), (0,1), (1,1).
#' d <- simulate_study(n = 500, prevalence = 0.2,
#'                     sensitivity = c(A = 0.9, B = 0.7),
#'                     specificity = c(A = 0.85, B = 0.95),
#'                     verification_prob = c(0.1, 0.9, 0.9, 1))
#' head(d)
#' # proportion of verified subjects by profile
#' tapply(d$verified, paste(d$A, d$B), mean)
#'
#' # Everyone verified: a single verification probability equal to 1
#' d_full <- simulate_study(n = 500, prevalence = 0.2,
#'                          sensitivity = c(A = 0.9, B = 0.7),
#'                          specificity = c(A = 0.85, B = 0.95),
#'                          verification_prob = 1)
#' mean(d_full$verified)
#'
#' # A never-verified profile: only subjects positive to at least one test are
#' # verified, so the disease status is missing for every (0,0) subject
#' d_screen <- simulate_study(n = 500, prevalence = 0.2,
#'                            sensitivity = c(A = 0.9, B = 0.7),
#'                            specificity = c(A = 0.85, B = 0.95),
#'                            verification_prob = c(0, 1, 1, 1))
#' table(paste(d_screen$A, d_screen$B), d_screen$verified)
#'
#' # Verification that depends on the disease (MNAR): with odds_ratio = 3 the
#' # odds of verification are three times higher for the diseased, so among
#' # the verified the diseased are over-represented
#' d_mnar <- simulate_study(n = 5000, prevalence = 0.2,
#'                          sensitivity = c(A = 0.9, B = 0.7),
#'                          specificity = c(A = 0.85, B = 0.95),
#'                          verification_prob = 0.3, odds_ratio = 3)
#' mean(d_mnar$disease, na.rm = TRUE)   # well above the prevalence 0.2
#'
#' # Three tests without names: the columns are called T1, T2, T3, and there
#' # are 2^3 = 8 profiles, so verification_prob has 8 values
#' d3 <- simulate_study(n = 500, prevalence = 0.1,
#'                      sensitivity = c(0.9, 0.8, 0.7),
#'                      specificity = c(0.9, 0.9, 0.95),
#'                      verification_prob = c(0.05, rep(0.8, 6), 1))
#' names(d3)
#'
#' # Three tests dependent given the disease: among the diseased tests A and B
#' # share a latent severity, C much less; among the healthy only A and B
#' # are (weakly) dependent
#' d_dep <- simulate_study(n = 5000, prevalence = 0.2,
#'                         sensitivity = c(A = 0.9, B = 0.8, C = 0.7),
#'                         specificity = c(A = 0.85, B = 0.9, C = 0.95),
#'                         verification_prob = c(0.05, rep(0.8, 6), 1),
#'                         dependence = latent_dependence(
#'                           diseased_loadings = c(1, 1, 0.3),
#'                           healthy_loadings = c(0.5, 0.5, 0)))
#' sick <- which(d_dep$disease == 1)
#' cor(d_dep[sick, c("A", "B", "C")])
#'
#' # Correlation 0.4 between every pair of tests among the diseased and -0.2
#' # among the healthy
#' d_probit <- simulate_study(n = 5000, prevalence = 0.2,
#'                            sensitivity = c(A = 0.9, B = 0.8, C = 0.7),
#'                            specificity = c(A = 0.85, B = 0.9, C = 0.95),
#'                            verification_prob = 1,
#'                            dependence = probit_dependence(
#'                              diseased_correlation = 0.4,
#'                              healthy_correlation = -0.2))
#' @seealso [independence()], [latent_dependence()], [probit_dependence()]
#' @export
simulate_study <- function(n, prevalence, sensitivity, specificity,
                           verification_prob, odds_ratio = 1,
                           dependence = independence()) {
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
  if (!S7::S7_inherits(dependence, test_dependence)) {
    stop("`dependence` must be a dependence mechanism, such as `independence()`.",
         call. = FALSE)
  }
  tests <- names(sensitivity)
  if (is.null(tests)) tests <- paste0("T", seq_len(K))
  if (any(c("verified", "disease") %in% tests)) {
    stop("Tests cannot be called `verified` or `disease`.", call. = FALSE)
  }

  disease <- stats::rbinom(n, 1L, prevalence)
  results <- draw_tests(dependence, disease = disease, sensitivity = sensitivity,
                        specificity = specificity)
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

# Dependence between the tests in simulate_study() -----------------------------

# Abstract parent of the dependence mechanisms; not exported.
test_dependence <- S7::new_class("test_dependence", abstract = TRUE)

# Draws the n x K matrix of 0/1 test results given the disease status.
draw_tests <- S7::new_generic(
  "draw_tests", "dependence",
  function(dependence, disease, sensitivity, specificity) S7::S7_dispatch()
)

#' Independent tests in simulated data
#'
#' `independence()` is the assumption, used when simulating data, that the
#' test results are independent given the disease status: among the diseased,
#' and separately among the healthy, the result of one test carries no
#' information on the results of the others. The function has no arguments.
#' The object is passed to [simulate_study()] as its `dependence` argument,
#' where it is the default.
#'
#' @details
#' Let \eqn{t = (t_1, \dots, t_K)} be a response profile, with \eqn{t_k = 1}
#' for a positive result, and let \eqn{Se_k} and \eqn{Sp_k} be the sensitivity
#' and specificity of test \eqn{k}. The probabilities of the profiles among the
#' diseased and among the healthy are the products
#' \deqn{P(T = t \mid D = 1) = \prod_{k=1}^{K} Se_k^{t_k} (1 - Se_k)^{1 - t_k},
#' \qquad P(T = t \mid D = 0) = \prod_{k=1}^{K} (1 - Sp_k)^{t_k} Sp_k^{1 - t_k}.}
#'
#' @return An object of class `independence`.
#' @seealso [simulate_study()], [latent_dependence()], [probit_dependence()]
#' @examples
#' independence()
#' d <- simulate_study(n = 100, prevalence = 0.2, sensitivity = c(0.9, 0.8),
#'                     specificity = c(0.9, 0.95), verification_prob = 1,
#'                     dependence = independence())
#' @export
independence <- S7::new_class("independence", parent = test_dependence)

S7::method(print_generic, independence) <- function(x, ...) {
  cat("<independence> tests independent given the disease status\n")
  invisible(x)
}

S7::method(draw_tests, independence) <- function(dependence, disease, sensitivity,
                                                 specificity) {
  # probability of a positive result: sensitivity if diseased, 1 - specificity if not
  p_positive <- outer(disease, sensitivity) + outer(1L - disease, 1 - specificity)
  matrix(stats::rbinom(length(p_positive), 1L, p_positive), length(disease))
}

#' Tests dependent through a shared latent variable in simulated data
#'
#' `latent_dependence()` makes the simulated test results dependent given the
#' disease status through one continuous latent variable shared by all the
#' tests of a subject, such as the severity of the condition among the
#' diseased or a condition that causes false positives among the healthy. The
#' sensitivity and specificity of each test do not change: the latent variable
#' changes only how often the tests are right or wrong together. The object is
#' passed to [simulate_study()] as its `dependence` argument.
#'
#' @details
#' Each subject has a latent variable \eqn{Z \sim N(0, 1)}, independent of the
#' disease status \eqn{D}. Given \eqn{D = d} and \eqn{Z = z} the \eqn{K} tests
#' are independent, with
#' \deqn{P(T_k = 1 \mid D = d, Z = z) = \Phi(a_{dk} + b_{dk} z),}
#' where \eqn{\Phi} is the standard normal distribution function and
#' \eqn{b_{dk}} is the loading of test \eqn{k} in class \eqn{d}
#' (`diseased_loadings` for \eqn{d = 1}, `healthy_loadings` for \eqn{d = 0}).
#' The intercepts are
#' \deqn{a_{1k} = \Phi^{-1}(Se_k)\sqrt{1 + b_{1k}^2}, \qquad
#' a_{0k} = \Phi^{-1}(1 - Sp_k)\sqrt{1 + b_{0k}^2},}
#' so that, integrating out \eqn{Z}, \eqn{P(T_k = 1 \mid D = 1) = Se_k} and
#' \eqn{P(T_k = 0 \mid D = 0) = Sp_k} exactly, for any value of the loadings.
#'
#' Equivalently, test \eqn{k} is positive when a normal latent variable
#' \eqn{Y_k = b_{dk} Z + \varepsilon_k}, with \eqn{\varepsilon_k \sim N(0, 1)}
#' independent, exceeds a threshold. Within class \eqn{d} the correlation of
#' \eqn{Y_j} and \eqn{Y_k} is
#' \deqn{\rho_{d,jk} = \frac{b_{dj}\, b_{dk}}{\sqrt{(1 + b_{dj}^2)(1 + b_{dk}^2)}}.}
#' Hence:
#' - a loading 0 makes its test independent of all the others in that class;
#' - tests \eqn{j} and \eqn{k} agree more often than under independence when
#'   \eqn{b_{dj} b_{dk} > 0}, and less often when \eqn{b_{dj} b_{dk} < 0};
#' - larger loadings in absolute value give stronger dependence. For example,
#'   two tests with sensitivity 0.8 and loadings 1 among the diseased have
#'   \eqn{\rho = 0.5} and an odds ratio of about 4.7 between their results;
#'   with loadings 2, \eqn{\rho = 0.8} and an odds ratio of about 19.
#'
#' **Allowed values.** Every combination of finite real values is allowed:
#' the loadings of a class are not constrained by each other, by the loadings
#' of the other class or by the sensitivities and specificities, and every
#' combination gives a valid distribution of the results with the requested
#' sensitivities and specificities. The two classes are independent choices,
#' so everything below holds separately among the diseased and among the
#' healthy.
#'
#' **Dependences that can be obtained.** Since all the tests share one latent
#' variable, the loadings cannot produce every pattern of pairwise
#' dependences. Write \eqn{s_k = |b_{dk}| / \sqrt{1 + b_{dk}^2}}, a number in
#' \eqn{[0, 1)} that measures how strongly test \eqn{k} depends on \eqn{Z};
#' then \eqn{|\rho_{d,jk}| = s_j s_k}. Hence:
#' - **strength:** each test has a single strength \eqn{s_k}, and the
#'   dependence of a pair is the product of the strengths of its two tests.
#'   For three tests \eqn{|\rho_{d,23}| = |\rho_{d,12}\, \rho_{d,13}| / s_1^2
#'   \ge |\rho_{d,12}\, \rho_{d,13}|}, so tests 2 and 3 cannot be much less
#'   dependent on each other than both are on test 1: for example
#'   \eqn{\rho_{d,12} = \rho_{d,13} = 0.8} forces \eqn{\rho_{d,23} \ge 0.64}.
#'   All the correlations are strictly between \eqn{-1} and 1;
#' - **sign:** the sign of \eqn{\rho_{d,jk}} is the product of the signs of
#'   \eqn{b_{dj}} and \eqn{b_{dk}}. The tests split into those with positive
#'   and those with negative loading, with positive dependence within each
#'   group and negative dependence between the groups. Changing the sign of all
#'   the loadings of a class gives the same distribution, so only the relative
#'   signs matter. For any three tests the product
#'   \eqn{\rho_{d,12}\, \rho_{d,13}\, \rho_{d,23}} is never negative: for
#'   example, three tests cannot all be negatively dependent on each other;
#' - **independence:** \eqn{\rho_{d,jk} = 0} only when \eqn{b_{dj} = 0} or
#'   \eqn{b_{dk} = 0}, and a loading 0 makes its test independent of all the
#'   others. Two tests cannot be independent of each other while both depend on
#'   a third.
#'
#' To give the same correlation \eqn{\rho \in [0, 1)} to every pair, use the
#' same loading \eqn{b = \sqrt{\rho / (1 - \rho)}} for all tests: \eqn{b = 0.5}
#' gives \eqn{\rho = 0.2}, \eqn{b = 1} gives \eqn{\rho = 0.5} and \eqn{b = 2}
#' gives \eqn{\rho = 0.8}. Patterns outside these limits, such as a correlation
#' of any sign and size for each pair, are obtained with
#' `dependence = probit_dependence(...)`.
#'
#' @param diseased_loadings Loadings \eqn{b_{1k}} of the tests among the
#'   diseased: a numeric vector of finite values, of any sign, either of length
#'   1 (the same loading for all tests) or of length \eqn{K} (one loading per
#'   test, in the order of the `sensitivity` argument of `simulate_study()`).
#'   Any combination of finite values is allowed; the Details describe the
#'   dependences each combination produces. Default `0`: the tests are
#'   independent among the diseased. The length is checked against \eqn{K}
#'   when the data are simulated.
#' @param healthy_loadings Loadings \eqn{b_{0k}} of the tests among the
#'   healthy, with the same form as `diseased_loadings`. Default `0`: the tests
#'   are independent among the healthy.
#'
#' @return An object of class `latent_dependence`.
#' @seealso [simulate_study()], [independence()], [probit_dependence()],
#'   [random_effects()], which fits the same latent structure
#' @examples
#' # same loading for all tests among the diseased, independence among the healthy
#' latent_dependence(diseased_loadings = 1)
#'
#' # one loading per test and class: A and B strongly dependent among the
#' # diseased, C almost independent; among the healthy only A and B, weakly
#' dep <- latent_dependence(diseased_loadings = c(1, 1, 0.3),
#'                          healthy_loadings = c(0.5, 0.5, 0))
#' d <- simulate_study(n = 2000, prevalence = 0.3,
#'                     sensitivity = c(A = 0.9, B = 0.8, C = 0.7),
#'                     specificity = c(A = 0.85, B = 0.9, C = 0.95),
#'                     verification_prob = 1, dependence = dep)
#' sick <- d$disease == 1
#' cor(d[sick, c("A", "B", "C")])
#'
#' # mixed signs: A and B agree, C disagrees with both
#' latent_dependence(diseased_loadings = c(1, 1, -1))
#' @export
latent_dependence <- S7::new_class(
  "latent_dependence",
  parent = test_dependence,
  properties = list(
    diseased_loadings = S7::class_numeric,
    healthy_loadings = S7::class_numeric
  ),
  constructor = function(diseased_loadings = 0, healthy_loadings = 0) {
    S7::new_object(S7::S7_object(), diseased_loadings = diseased_loadings,
                   healthy_loadings = healthy_loadings)
  },
  validator = function(self) {
    for (name in c("diseased_loadings", "healthy_loadings")) {
      x <- S7::prop(self, name)
      if (length(x) == 0L || !all(is.finite(x))) {
        return(paste0("@", name, " must be one or more finite numbers"))
      }
    }
  }
)

S7::method(print_generic, latent_dependence) <- function(x, ...) {
  cat("<latent_dependence> tests dependent through a shared latent variable\n",
      "  diseased_loadings: ", paste(format(x@diseased_loadings), collapse = ", "), "\n",
      "  healthy_loadings: ", paste(format(x@healthy_loadings), collapse = ", "), "\n",
      sep = "")
  invisible(x)
}

# Recycles a per-test value of length 1 to the K tests, or checks its length.
per_test <- function(x, name, K) {
  if (!length(x) %in% c(1L, K)) {
    stop("`", name, "` must have length 1 or one value per test (", K, ").",
         call. = FALSE)
  }
  rep_len(x, K)
}

S7::method(draw_tests, latent_dependence) <- function(dependence, disease, sensitivity,
                                                      specificity) {
  K <- length(sensitivity)
  b <- rbind(per_test(dependence@healthy_loadings, "healthy_loadings", K),
             per_test(dependence@diseased_loadings, "diseased_loadings", K))
  # marginal probits by class; row 1 healthy, row 2 diseased (indexed by disease + 1)
  m <- rbind(stats::qnorm(1 - specificity), stats::qnorm(sensitivity))
  b <- b[disease + 1L, , drop = FALSE]
  a <- m[disease + 1L, , drop = FALSE] * sqrt(1 + b^2)
  z <- stats::rnorm(length(disease))
  p_positive <- stats::pnorm(a + b * z)
  matrix(stats::rbinom(length(p_positive), 1L, p_positive), length(disease))
}

#' Tests with correlated latent variables in simulated data
#'
#' `probit_dependence()` makes the simulated test results dependent given the
#' disease status by giving each test a normal latent variable and letting
#' these variables be correlated, with one correlation for each pair of tests
#' and disease class. The correlations can have any sign, separately for each
#' pair. The sensitivity and specificity of each test do not change: the
#' correlations change only how often the tests are right or wrong together.
#' The object is passed to [simulate_study()] as its `dependence` argument.
#'
#' @details
#' Given the disease status \eqn{D = d}, the \eqn{K} tests have latent
#' variables \eqn{Y = (Y_1, \dots, Y_K)} jointly normal with mean 0, variance 1
#' and correlation matrix \eqn{R_d} (`diseased_correlation` for \eqn{d = 1},
#' `healthy_correlation` for \eqn{d = 0}). Test \eqn{k} is positive when
#' \deqn{Y_k < \Phi^{-1}(Se_k) \textrm{ if } d = 1, \qquad
#' Y_k < \Phi^{-1}(1 - Sp_k) \textrm{ if } d = 0,}
#' where \eqn{\Phi} is the standard normal distribution function, so that
#' \eqn{P(T_k = 1 \mid D = 1) = Se_k} and \eqn{P(T_k = 0 \mid D = 0) = Sp_k}
#' exactly, whatever the correlations.
#'
#' The element \eqn{R_{d,jk}} is the correlation between the latent variables
#' of tests \eqn{j} and \eqn{k} in class \eqn{d} (known as the tetrachoric
#' correlation of the two results). Positive values make the two tests agree
#' more often than under independence, negative values less often, and 0 makes
#' their results uncorrelated. For example, two tests with sensitivity 0.8 and
#' correlation 0.5 among the diseased have an odds ratio of about 4.7 between
#' their results; with correlation 0.8, about 19. With all the correlations of
#' a class equal to 0 the tests are independent in that class.
#'
#' A correlation matrix must be positive definite: not every set of pairwise
#' correlations is possible. For example, the same correlation \eqn{r} for all
#' pairs of \eqn{K} tests needs \eqn{-1/(K - 1) < r < 1}, so with three tests
#' \eqn{r > -1/2}. Since \eqn{K} is known only when the data are simulated,
#' this condition is checked by `simulate_study()`, which gives an error for a
#' matrix that is not positive definite.
#'
#' @param diseased_correlation Correlations \eqn{R_1} between the latent
#'   variables of the tests among the diseased: either a single number, the
#'   same correlation for every pair of tests, or a symmetric
#'   \eqn{K \times K} matrix with 1 on the diagonal, whose rows and columns
#'   follow the order of the `sensitivity` argument of `simulate_study()`.
#'   Values in \eqn{[-1, 1]}. Default `0`: the tests are independent among the
#'   diseased. The size of the matrix and its positive definiteness (which
#'   excludes, for example, a single number equal to 1 with two or more tests)
#'   are checked when the data are simulated.
#' @param healthy_correlation Correlations \eqn{R_0} between the latent
#'   variables of the tests among the healthy, with the same form as
#'   `diseased_correlation`. Default `0`: the tests are independent among the
#'   healthy.
#'
#' @return An object of class `probit_dependence`.
#' @seealso [simulate_study()], [independence()], [latent_dependence()]
#' @examples
#' # correlation 0.4 between every pair among the diseased, 0.2 among the healthy
#' probit_dependence(diseased_correlation = 0.4, healthy_correlation = 0.2)
#'
#' # different correlation for each pair: A and B positively correlated,
#' # C negatively correlated with both
#' R1 <- matrix(c( 1,    0.5, -0.3,
#'                 0.5,  1,   -0.3,
#'                -0.3, -0.3,  1), 3, 3)
#' d <- simulate_study(n = 2000, prevalence = 0.3,
#'                     sensitivity = c(A = 0.9, B = 0.8, C = 0.7),
#'                     specificity = c(A = 0.85, B = 0.9, C = 0.95),
#'                     verification_prob = 1,
#'                     dependence = probit_dependence(diseased_correlation = R1))
#' sick <- d$disease == 1
#' cor(d[sick, c("A", "B", "C")])
#'
#' # every pair negatively correlated: possible with three tests if r > -1/2
#' probit_dependence(diseased_correlation = -0.3)
#' @export
probit_dependence <- S7::new_class(
  "probit_dependence",
  parent = test_dependence,
  properties = list(
    diseased_correlation = S7::class_numeric,
    healthy_correlation = S7::class_numeric
  ),
  constructor = function(diseased_correlation = 0, healthy_correlation = 0) {
    S7::new_object(S7::S7_object(), diseased_correlation = diseased_correlation,
                   healthy_correlation = healthy_correlation)
  },
  validator = function(self) {
    for (name in c("diseased_correlation", "healthy_correlation")) {
      x <- S7::prop(self, name)
      if (length(x) == 0L || anyNA(x) || any(abs(x) > 1)) {
        return(paste0("@", name, " must have values in [-1, 1]"))
      }
      if (is.null(dim(x))) {
        if (length(x) != 1L) {
          return(paste0("@", name, " must be a single number or a correlation matrix"))
        }
      } else if (!is.matrix(x) || nrow(x) != ncol(x) || !isSymmetric(unname(x)) ||
                 any(diag(x) != 1)) {
        return(paste0("@", name, " must be a symmetric matrix with 1 on the diagonal"))
      }
    }
  }
)

S7::method(print_generic, probit_dependence) <- function(x, ...) {
  label <- function(r) {
    if (is.matrix(r)) paste0(nrow(r), " x ", ncol(r), " matrix") else format(r)
  }
  cat("<probit_dependence> tests with correlated latent variables\n",
      "  diseased_correlation: ", label(x@diseased_correlation), "\n",
      "  healthy_correlation: ", label(x@healthy_correlation), "\n", sep = "")
  invisible(x)
}

# Cholesky factor of the correlation matrix for K tests, given a single number
# or a K x K matrix; stops if it does not fit K or is not positive definite.
correlation_factor <- function(r, name, K) {
  if (!is.matrix(r)) {
    r <- matrix(r, K, K)
    diag(r) <- 1
  } else if (nrow(r) != K) {
    stop("`", name, "` must be a single number or a ", K, " x ", K, " matrix.",
         call. = FALSE)
  }
  tryCatch(chol(r), error = function(e) {
    stop("`", name, "` is not a valid correlation matrix for ", K,
         " tests (it is not positive definite).", call. = FALSE)
  })
}

S7::method(draw_tests, probit_dependence) <- function(dependence, disease, sensitivity,
                                                      specificity) {
  K <- length(sensitivity)
  factors <- list(correlation_factor(dependence@healthy_correlation,
                                     "healthy_correlation", K),
                  correlation_factor(dependence@diseased_correlation,
                                     "diseased_correlation", K))
  # a test is positive when its latent variable is below this threshold;
  # row 1 healthy, row 2 diseased
  thresholds <- rbind(stats::qnorm(1 - specificity), stats::qnorm(sensitivity))
  results <- matrix(0L, length(disease), K)
  for (d in 0:1) {
    rows <- which(disease == d)
    y <- matrix(stats::rnorm(length(rows) * K), length(rows), K) %*% factors[[d + 1L]]
    results[rows, ] <- 1L * (y < rep(thresholds[d + 1L, ], each = length(rows)))
  }
  results
}
