#' @include priors.R verification.R generics.R
NULL

# Accuracy models ----------------------------------------------------------------
#
# A model object holds the priors and the verification mechanism. The checks
# that depend on the data (for example, never-verified profiles that need a
# prior) are done by baysesp().

#' Accuracy models
#'
#' `accuracy_model` is the abstract parent class of the models for the joint
#' distribution of the test results and of the disease status. It cannot be
#' created directly: the concrete models are [saturated()], [loglinear()] and
#' [random_effects()], each with its own constructor.
#'
#' Every model has the property `verification`, an object of class
#' `verification_mechanism` (`mar()` or `mnar()`) that states how the
#' probability of being verified depends on the test results and on the
#' disease status. A model object holds the priors and the verification
#' mechanism; the data are supplied to `baysesp()`, which fits the model.
#'
#' @seealso [saturated()], [loglinear()], [random_effects()]
#' @export
accuracy_model <- S7::new_class(
  "accuracy_model",
  abstract = TRUE,
  properties = list(verification = verification_mechanism)
)

#' Saturated model
#'
#' The model with no assumption on the dependence between the tests: the
#' distribution of the profiles is free within each disease class. It is
#' written in a transparent parametrization, which separates the identified
#' parameters from the non-identified ones:
#' - identified: the frequency of each response profile, the verification
#'   probability of each profile with at least one verified subject, and the
#'   proportion of diseased among its verified subjects;
#' - not identified: the proportion of diseased among the unverified subjects
#'   of each such profile, which is fixed by the verification mechanism, and
#'   the proportion of diseased in each never-verified profile, whose posterior
#'   equals its prior (Gustafson, 2015).
#'
#' @details
#' With \eqn{K} tests there are \eqn{J = 2^K} response profiles \eqn{t}. For
#' each profile the model has four parameters:
#' - \eqn{\omega_t = P(T = t)}: frequency of the profile among all subjects;
#' - \eqn{v_t = P(V = 1 \mid T = t)}: probability that a subject with profile
#'   \eqn{t} is verified;
#' - \eqn{\phi^{(1)}_t = P(D = 1 \mid T = t, V = 1)}: proportion of diseased
#'   among the verified subjects with profile \eqn{t};
#' - \eqn{\phi^{(0)}_t = P(D = 1 \mid T = t, V = 0)}: proportion of diseased
#'   among the unverified subjects with profile \eqn{t}.
#'
#' Let \eqn{n_{t,1}} and \eqn{n_{t,0}} be the numbers of verified diseased and
#' verified healthy subjects with profile \eqn{t}, \eqn{u_t} the number of
#' unverified subjects, and \eqn{N_t = n_{t,1} + n_{t,0} + u_t}. The
#' log-likelihood is
#' \deqn{\ell = \sum_t \Big[N_t \log \omega_t + (n_{t,1} + n_{t,0}) \log v_t + u_t \log(1 - v_t) + n_{t,1} \log \phi^{(1)}_t + n_{t,0} \log(1 - \phi^{(1)}_t)\Big].}
#' The unverified subjects enter only through \eqn{N_t} and \eqn{u_t}, because
#' their disease status is not observed (Begg and Greenes, 1983; Zhou, 1993). The proportion \eqn{\phi^{(0)}_t}
#' does not enter the log-likelihood: the data carry no information on it (Gustafson, 2015), so
#' it is determined by the verification mechanism and by \eqn{\phi^{(1)}_t}.
#' A never-verified profile has no verified subject, so \eqn{v_t = 0} and only
#' its first term remains.
#'
#' With `verification = mar()` (missing at random) verification depends only on the test
#' results: among subjects with the same profile, verified and unverified
#' subjects have the same proportion of diseased, so
#' \eqn{\phi^{(0)}_t = \phi^{(1)}_t}.
#'
#' With `verification = mnar()` (missing not at random) verification may also depend on the
#' disease status (Kosinski and Barnhart, 2003a):
#' \deqn{\mathrm{logit}\,\phi^{(0)}_t = \mathrm{logit}\,\phi^{(1)}_t - \beta,}
#' where \eqn{\beta} is the log odds ratio of verification between diseased and
#' healthy subjects with the same profile. Equivalently, \eqn{\exp(-\beta)} is
#' the odds ratio of disease between unverified and verified subjects with the
#' same profile: \eqn{\beta > 0} means that the unverified are less often
#' diseased. There is one \eqn{\beta} common to all profiles with verified
#' subjects, or one per profile with `verification = mnar(by_profile = TRUE)`. Since
#' \eqn{\beta} does not enter the log-likelihood, its posterior equals its
#' prior, and the results depend on the assumed value: they should be examined
#' over a range of plausible values (Gustafson, 2015; Kosinski and Barnhart, 2003b).
#'
#' The quantities of interest are obtained as follows. For each profile,
#' \deqn{\psi_t = v_t \phi^{(1)}_t + (1 - v_t) \phi^{(0)}_t}
#' is the proportion of diseased among subjects with profile \eqn{t}; for a
#' never-verified profile \eqn{\phi^{(0)}_t = q_t} is the proportion of diseased
#' in it, so \eqn{\psi_t = q_t}. Then
#' \deqn{\pi = \sum_t \omega_t \psi_t, \qquad \theta_{1,t} = \frac{\omega_t \psi_t}{\pi}, \qquad \theta_{0,t} = \frac{\omega_t (1 - \psi_t)}{1 - \pi},}
#' where \eqn{\pi} is the prevalence, \eqn{\theta_{1,t} = P(T = t \mid D = 1)}
#' and \eqn{\theta_{0,t} = P(T = t \mid D = 0)}. The sensitivity of a test is
#' the sum of \eqn{\theta_{1,t}} over the profiles with a positive result for
#' that test, and its specificity the sum of \eqn{\theta_{0,t}} over the
#' profiles with a negative result. Under both mechanisms the never-verified
#' profiles enter \eqn{\pi} through \eqn{q_t}, whose posterior equals its prior.
#'
#' @param frequency_prior [prior_dirichlet()] on the vector of profile
#'   frequencies \eqn{\omega = (\omega_1, \dots, \omega_J)}. Its
#'   `concentration` is a single positive number, the same for every profile,
#'   or a vector with one positive value per profile, in binary order: the
#'   profile with all tests negative comes first, and the first test changes
#'   fastest (with two tests: (0,0), (1,0), (0,1), (1,1)). Default
#'   `prior_dirichlet(1)`, the uniform distribution on the simplex.
#' @param verified_prior Prior on the verification probability \eqn{v_t} of each
#'   profile with at least one verified subject. The same prior is used,
#'   independently, for each such profile. A prior on a probability, such as
#'   [prior_beta()] or [elicit_beta()], acts on \eqn{v_t} itself; a prior on the
#'   real line, such as [prior_normal()] or [prior_student()], acts on
#'   \eqn{\mathrm{logit}\,v_t}. Default `prior_beta(0.5, 2)`, the uniform
#'   distribution on (0, 1). Never-verified profiles have no such parameter.
#' @param profile_prevalence_prior Prior on \eqn{\phi^{(1)}_t}, the proportion of
#'   diseased among the verified subjects of each profile with at least one
#'   verified subject. The same prior is used, independently, for each such
#'   profile, with the same rules as `verified_prior`: a prior on a probability
#'   acts on the proportion, a prior on the real line on its logit. Default
#'   `prior_beta(0.5, 1)`, the Jeffreys prior Beta(1/2, 1/2). It matters in
#'   profiles with few verified subjects and almost all of one class, typically
#'   the profile with all tests negative when few negatives are verified: there
#'   the uniform `prior_beta(0.5, 2)` pulls \eqn{\phi^{(1)}_t} towards 0.5 more,
#'   moving diseased subjects into that profile, so that the sensitivities are
#'   underestimated. In simulations with 5% to 10% of the negatives verified,
#'   the Jeffreys prior about halved this bias. With this prior the posterior of
#'   \eqn{\mathrm{logit}\,\phi^{(1)}_t} can have a long tail, so `baysesp()`
#'   uses `adapt_delta = 0.95` for this model.
#' @param unverified_prevalence_prior Prior on \eqn{q_t}, the proportion of
#'   diseased in each never-verified profile (a profile with no verified
#'   subject). It has no default prior: its default `NULL` means that no prior
#'   is given. It is required when at least one profile is never verified, and
#'   `baysesp()` stops with an error if it is missing in that case. Its value
#'   can be one prior, which applies to every never-verified profile, or a list
#'   with one prior per never-verified profile, in binary order among those
#'   profiles (with two tests, if both (0,1) and (1,1) are never verified, the
#'   list has two elements: the first for (0,1), the second for (1,1)). A list
#'   of length one is the same as a single prior. Each prior has the same
#'   allowed classes as `verified_prior`. Its posterior equals its prior, so it
#'   is an explicit assumption; use `elicit_beta()` to specify it from a
#'   plausible value and an upper bound.
#' @param verification Verification mechanism: `mar()` (the default), where
#'   verification depends only on the test results, or `mnar()`, where it also
#'   depends on the disease status. See Details for what changes.
#' @return An object of class `saturated`.
#' @references
#' Begg, C. B. and Greenes, R. A. (1983). Assessment of diagnostic tests when disease verification is subject to selection bias. *Biometrics*, 39(1), 207–215. \doi{10.2307/2530820}
#'
#' Gustafson, P. (2015). *Bayesian Inference for Partially Identified Models: Exploring the Limits of Limited Data*. Chapman & Hall/CRC, Boca Raton, FL.
#'
#' Kosinski, A. S. and Barnhart, H. X. (2003a). Accounting for nonignorable verification bias in assessment of diagnostic tests. *Biometrics*, 59(1), 163–171. \doi{10.1111/1541-0420.00019}
#'
#' Kosinski, A. S. and Barnhart, H. X. (2003b). A global sensitivity analysis of performance of a medical diagnostic test when verification bias is present. *Statistics in Medicine*, 22(17), 2711–2721. \doi{10.1002/sim.1517}
#'
#' Zhou, X.-H. (1993). Maximum likelihood estimators of sensitivity and specificity corrected for verification bias. *Communications in Statistics -- Theory and Methods*, 22(11), 3177–3198. \doi{10.1080/03610929308831209}
#' @examples
#' saturated()
#' saturated(unverified_prevalence_prior = elicit_beta(mean = 0.005, upper = 0.015))
#' saturated(verification = mnar())
#' @export
saturated <- S7::new_class(
  "saturated",
  parent = accuracy_model,
  properties = list(
    frequency_prior = prior_dirichlet,
    verified_prior = prior_distribution,
    profile_prevalence_prior = prior_distribution,
    unverified_prevalence_prior = S7::new_union(NULL, prior_distribution, S7::class_list)
  ),
  constructor = function(frequency_prior = prior_dirichlet(1),
                         verified_prior = prior_beta(0.5, 2),
                         profile_prevalence_prior = prior_beta(0.5, 1),
                         unverified_prevalence_prior = NULL,
                         verification = mar()) {
    S7::new_object(
      S7::S7_object(),
      verification = verification,
      frequency_prior = frequency_prior,
      verified_prior = verified_prior,
      profile_prevalence_prior = profile_prevalence_prior,
      unverified_prevalence_prior = unverified_prevalence_prior
    )
  },
  validator = function(self) {
    priors <- c(list(verified_prior = self@verified_prior,
                     profile_prevalence_prior = self@profile_prevalence_prior),
                as_prior_list(self@unverified_prevalence_prior,
                              "unverified_prevalence_prior"))
    for (name in names(priors)) {
      problem <- tryCatch(check_prior(priors[[name]], "probability", name),
                          error = conditionMessage)
      if (is.character(problem)) return(problem)
    }
  }
)

# A prior or a list of priors as a named list (names for error messages).
as_prior_list <- function(x, argument) {
  if (is.null(x)) return(list())
  if (!is.list(x)) x <- list(x)
  stats::setNames(x, paste0(argument, if (length(x) > 1L) paste0("[[", seq_along(x), "]]")))
}

S7::method(print_generic, saturated) <- function(x, ...) {
  q <- x@unverified_prevalence_prior
  q <- if (is.null(q)) "none" else
    paste(vapply(as_prior_list(q, ""), prior_label, character(1)), collapse = "; ")
  cat("<saturated> saturated model, transparent parametrization\n",
      "  frequency_prior:             ", prior_label(x@frequency_prior), "\n",
      "  verified_prior:              ", prior_label(x@verified_prior), "\n",
      "  profile_prevalence_prior:    ", prior_label(x@profile_prevalence_prior), "\n",
      "  unverified_prevalence_prior: ", q, "\n",
      "  verification:                ", sub("^baysesp::", "", class(x@verification)[1]),
      "\n", sep = "")
  invisible(x)
}
