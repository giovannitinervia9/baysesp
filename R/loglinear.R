#' @include models.R generics.R
NULL

#' Log-linear model of the test profiles
#'
#' Defines a model in which the distribution of the response profiles is
#' log-linear in each disease class, so the dependence between the tests is
#' restricted to interactions of at most `interaction_order` tests. The
#' object is passed to [baysesp()] as its `model` argument.
#'
#' @details
#' Let \eqn{K} be the number of tests, \eqn{t = (t_1, \dots, t_K)} a response
#' profile with \eqn{t_k = 1} for a positive result, \eqn{d = 1} for diseased
#' and \eqn{d = 0} for healthy, and let \eqn{m} be `interaction_order`. For each
#' disease class the profile probabilities \eqn{\theta_d} are
#' \deqn{\theta_d = \mathrm{softmax}(X_m \eta_d),}
#' where \eqn{X_m} is the design matrix whose columns are the functions
#' \eqn{x_S(t) = \prod_{k \in S} t_k}, one for each set \eqn{S} of at most
#' \eqn{m} tests, and \eqn{\eta_d} is the vector of the corresponding
#' class-specific parameters \eqn{\eta_{d,S}}. Equivalently,
#' \deqn{\theta_{d,t} \propto \exp\Big(\sum_{S} \eta_{d,S}\, x_S(t)\Big),}
#' where the normalisation runs over all \eqn{2^K} profiles.
#'
#' The order controls the structure of the dependence:
#' \itemize{
#'   \item Order 1 (`interaction_order = 1`): the tests are conditionally
#'     independent given the disease status.
#'   \item Order 2: each pair of tests may be associated within each class,
#'     but there are no interactions among three or more tests.
#'   \item Order \eqn{K}: no restriction, every profile distribution is
#'     allowed in each class.
#' }
#'
#' The parameters have a direct interpretation. With order 1,
#' \eqn{\eta_{1,k} = \mathrm{logit}(Se_k)} and
#' \eqn{\eta_{0,k} = \mathrm{logit}(1 - Sp_k)}, so the main effects are the
#' logits of the sensitivity and of the false positive rate of test \eqn{k}.
#' A pairwise term \eqn{\eta_{d,\{j,k\}}} is the log odds ratio between tests
#' \eqn{j} and \eqn{k} within class \eqn{d}, given the results of the other
#' tests (with order 2 it is the same whatever those results; with higher
#' orders it is the value when the other tests are negative): zero means that the two tests are
#' independent within that class, positive values mean that they are positively
#' associated among the diseased (\eqn{d = 1}) or among the healthy (\eqn{d = 0}).
#'
#' The log-likelihood of the data, with \eqn{n_{t,d}} verified subjects of
#' profile \eqn{t} and disease status \eqn{d}, \eqn{u_t} unverified subjects of
#' profile \eqn{t}, \eqn{\pi_1 = \pi} the prevalence and \eqn{\pi_0 = 1 - \pi},
#' is
#' \deqn{\ell = \sum_t \Big[\sum_d n_{t,d} \log(\pi_d\,\theta_{d,t}\,\lambda_{t,d}) + u_t \log \sum_d \pi_d\,\theta_{d,t}\,(1 - \lambda_{t,d})\Big],}
#' where \eqn{\lambda_{t,d}} is the probability of verification.
#' With `verification = mar()` the verification probability does not depend on the disease
#' status given the profile, so the verification term is separated from the
#' disease model and can be ignored:
#' \deqn{\ell = \sum_t \Big[\sum_d n_{t,d} \log(\pi_d\,\theta_{d,t}) + u_t \log \sum_d \pi_d\,\theta_{d,t}\Big]}
#' plus a term that depends only on the verification probabilities.
#' With `verification = mnar()` the verification probability may depend on the disease status:
#' \deqn{\mathrm{logit}\,\lambda_{t,d} = \alpha_t + \beta d,}
#' with one \eqn{\alpha_t} for each verified profile. The parameter \eqn{\beta}
#' is the log odds ratio of verification between diseased and healthy subjects
#' with the same profile, and \eqn{\exp(-\beta)} is the odds ratio of disease
#' between unverified and verified subjects with the same profile. With
#' `verification = mnar(by_profile = TRUE)` there is one \eqn{\beta_t} for each profile.
#' Profiles that are never verified have \eqn{\lambda_{t,d} = 0}.
#'
#' With order lower than \eqn{K} the restricted structure also determines the
#' disease proportion in the never-verified profiles, and makes \eqn{\beta}
#' formally estimable. This information comes from the assumed structure, not
#' from the data, and it is fragile if the structure is wrong. With order
#' \eqn{K} the model imposes no restriction and \eqn{\beta} is not identified.
#'
#' @param interaction_order A single whole number from 1 to the number of tests
#'   \eqn{K}, the largest number of tests in an interaction. Default `1L`
#'   (conditional independence). The lower bound is checked when the model is
#'   created; the upper bound is checked when the model is fitted.
#' @param prevalence_prior Prior on the prevalence \eqn{\pi}: a single prior,
#'   either on a probability, such as `prior_beta()`, which acts on \eqn{\pi},
#'   or on the real line, such as `prior_normal()`, which acts on
#'   \eqn{\mathrm{logit}\,\pi}. Default `prior_beta(0.5, 2)`, the uniform
#'   distribution on (0, 1).
#' @param main_effects_prior A single prior on the real line, applied to each
#'   class-specific main effect \eqn{\eta_{d,k}} (one per test and disease
#'   class). Default `prior_normal(0, 1.5)`. With order 1 these parameters are
#'   logits of sensitivities and of false positive rates.
#' @param interactions_prior A single prior on the real line, applied to each
#'   class-specific interaction of order 2 or higher. Default
#'   `prior_normal(0, 1)`. It has no effect when `interaction_order = 1`.
#' @param verification The verification mechanism: `mar()` (default), with
#'   verification that depends only on the test profile, or `mnar()`, with
#'   verification that may also depend on the disease status.
#' @return An object of class `loglinear`.
#' @seealso [saturated()] for the model without restrictions on the profile
#'   distribution, [mar()] and [mnar()] for the verification mechanisms.
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

S7::method(print_generic, loglinear) <- function(x, ...) {
  cat("<loglinear> log-linear model of order ", x@interaction_order, "\n",
      "  prevalence_prior:   ", prior_label(x@prevalence_prior), "\n",
      "  main_effects_prior: ", prior_label(x@main_effects_prior), "\n",
      "  interactions_prior: ", prior_label(x@interactions_prior), "\n",
      "  verification:       ", sub("^baysesp::", "", class(x@verification)[1]),
      "\n", sep = "")
  invisible(x)
}
