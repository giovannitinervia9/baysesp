#' @include priors.R verification.R
NULL

# Accuracy models ----------------------------------------------------------------
#
# A model object holds the priors and the verification mechanism. The checks
# that depend on the data (for example, never-verified profiles that need a
# prior) are done by baysesp().

#' Accuracy models
#'
#' `accuracy_model` is the abstract parent class of the models for the joint
#' distribution of the tests. Every model has a `verification` property, the
#' assumed [verification_mechanism].
#'
#' @param verification Verification mechanism, by default [mar()].
#' @export
accuracy_model <- S7::new_class(
  "accuracy_model",
  abstract = TRUE,
  properties = list(verification = verification_mechanism)
)

#' Saturated model
#'
#' The model of the main analysis (road 1 in the specification): no assumption
#' on the dependence between the tests. It is written in the transparent
#' parametrization, which separates the parameters the data identify from the
#' assumptions they cannot check:
#' - identified: the frequency of each response profile, the fraction of
#'   verified subjects in each profile and the proportion of diseased among
#'   them;
#' - not identified: the proportion of diseased in the profiles that are never
#'   verified. Its posterior is equal to its prior, so the prior is an explicit
#'   assumption and has no default.
#'
#' Priors with support `"real"`, such as [prior_normal()], act on the logit of
#' the probability they are given to.
#'
#' @param frequency_prior [prior_dirichlet()] on the frequencies of the
#'   response profiles.
#' @param verified_prior Prior on the fraction of verified subjects in each
#'   profile.
#' @param profile_prevalence_prior Prior on the proportion of diseased among
#'   the verified subjects of each profile.
#' @param unverified_prevalence_prior Prior on the proportion of diseased in
#'   each never-verified profile: one prior for all of them, or a list with one
#'   prior per never-verified profile, in binary order. Required when some
#'   profile has no verified subject; see [elicit_beta()].
#' @param verification Verification mechanism, by default [mar()].
#' @return An object of class `saturated`.
#' @examples
#' saturated()
#' saturated(unverified_prevalence_prior = elicit_beta(mean = 0.005, upper = 0.015))
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
                         profile_prevalence_prior = prior_beta(0.5, 2),
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

S7::method(print, saturated) <- function(x, ...) {
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
