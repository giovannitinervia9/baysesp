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

#' Log-linear models
#'
#' Models for comparison (road 2 in the specification): they assume a
#' structure for the dependence between the tests. The distribution of the
#' profiles in each disease class is log-linear with the parameters of the
#' subsets of at most `interaction_order` tests: order 1 is conditional
#' independence, order 2 allows the association of each pair of tests, and
#' order `K` (the number of tests) is the saturated model.
#'
#' With order lower than `K`, the structure gives information also on the
#' proportion of diseased in never-verified profiles, and makes the odds ratio
#' of [mnar()] verification formally estimable. That information comes from
#' the assumption, not from the data: use these models as a comparison with
#' [saturated()].
#'
#' @param interaction_order Highest number of tests in an interaction.
#' @param prevalence_prior Prior on the prevalence.
#' @param main_effects_prior Prior on the main effects (one per test and class).
#'   With order 1 the main effects are the logits of the sensitivities and of
#'   one minus the specificities.
#' @param interactions_prior Prior on the interactions.
#' @param verification Verification mechanism, by default [mar()].
#' @return An object of class `loglinear`.
#' @examples
#' loglinear()                       # conditional independence
#' loglinear(interaction_order = 2)  # pairwise dependence
#' @export
loglinear <- S7::new_class(
  "loglinear",
  parent = accuracy_model,
  properties = list(
    interaction_order = S7::class_integer,
    prevalence_prior = prior_distribution,
    main_effects_prior = prior_distribution,
    interactions_prior = prior_distribution
  ),
  constructor = function(interaction_order = 1L,
                         prevalence_prior = prior_beta(0.5, 2),
                         main_effects_prior = prior_normal(0, 1.5),
                         interactions_prior = prior_normal(0, 1),
                         verification = mar()) {
    if (!is_scalar_number(interaction_order) || interaction_order < 1 ||
        interaction_order != round(interaction_order)) {
      stop("`interaction_order` must be a positive whole number.", call. = FALSE)
    }
    S7::new_object(
      S7::S7_object(),
      verification = verification,
      interaction_order = as.integer(interaction_order),
      prevalence_prior = prevalence_prior,
      main_effects_prior = main_effects_prior,
      interactions_prior = interactions_prior
    )
  },
  validator = function(self) {
    targets <- c(prevalence_prior = "probability", main_effects_prior = "real",
                 interactions_prior = "real")
    for (name in names(targets)) {
      problem <- tryCatch(check_prior(S7::prop(self, name), targets[[name]], name),
                          error = conditionMessage)
      if (is.character(problem)) return(problem)
    }
  }
)

S7::method(print, loglinear) <- function(x, ...) {
  cat("<loglinear> log-linear model of order ", x@interaction_order, "\n",
      "  prevalence_prior:   ", prior_label(x@prevalence_prior), "\n",
      "  main_effects_prior: ", prior_label(x@main_effects_prior), "\n",
      "  interactions_prior: ", prior_label(x@interactions_prior), "\n",
      "  verification:       ", sub("^baysesp::", "", class(x@verification)[1]),
      "\n", sep = "")
  invisible(x)
}
