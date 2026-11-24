#' @include priors.R
NULL

#' Verification mechanisms
#'
#' `verification_mechanism` is the abstract parent class of the assumptions on
#' how subjects are selected for verification by the reference standard.
#' `mar()` assumes that verification depends only on the test results
#' (missing at random): among subjects with the same profile, verified and
#' unverified subjects have the same proportion of diseased.
#'
#' @return An object of class `mar`.
#' @examples
#' mar()
#' @export
verification_mechanism <- S7::new_class("verification_mechanism", abstract = TRUE)

#' @rdname verification_mechanism
#' @export
mar <- S7::new_class("mar", parent = verification_mechanism)

S7::method(print, mar) <- function(x, ...) {
  cat("<mar> verification depends only on the test results\n")
  invisible(x)
}

#' Verification that depends on the disease
#'
#' `mnar()` assumes that verification depends on the test results and also on
#' the disease status (missing not at random):
#' \eqn{logit P(V = 1 | T = t, D = d) = \alpha_t + \beta d}. The parameter
#' \eqn{\beta} is the log odds ratio of verification between diseased and
#' healthy subjects with the same profile; equivalently, \eqn{\exp(-\beta)} is
#' the odds ratio of disease between unverified and verified subjects with the
#' same profile (Proposition 7.2 of the specification). With \eqn{\beta > 0}
#' the diseased are verified more often, and the unverified are less often
#' diseased than MAR would say.
#'
#' In the saturated model \eqn{\beta} does not enter the likelihood: its
#' posterior equals its prior, which is therefore an explicit assumption. Use
#' [sensitivity_analysis()] to see how the results change over a range of
#' odds ratios.
#'
#' @param odds_ratio_prior Prior on \eqn{\beta}, the log odds ratio of
#'   verification between diseased and healthy subjects.
#' @param intercepts_prior Prior on the intercepts \eqn{\alpha_t}, used only by
#'   the structured models.
#' @param by_profile If `TRUE`, a different \eqn{\beta} for each profile, all
#'   with the same prior.
#' @return An object of class `mnar`.
#' @examples
#' mnar()
#' mnar(odds_ratio_prior = prior_normal(mean = log(2), sd = 0.5))
#' @export
mnar <- S7::new_class(
  "mnar",
  parent = verification_mechanism,
  properties = list(
    odds_ratio_prior = prior_distribution,
    intercepts_prior = prior_distribution,
    by_profile = S7::class_logical
  ),
  constructor = function(odds_ratio_prior = prior_normal(0, 1),
                         intercepts_prior = prior_normal(0, 2.5),
                         by_profile = FALSE) {
    S7::new_object(S7::S7_object(), odds_ratio_prior = odds_ratio_prior,
                   intercepts_prior = intercepts_prior, by_profile = by_profile)
  },
  validator = function(self) {
    for (name in c("odds_ratio_prior", "intercepts_prior")) {
      problem <- tryCatch(check_prior(S7::prop(self, name), "real", name),
                          error = conditionMessage)
      if (is.character(problem)) return(problem)
    }
    if (length(self@by_profile) != 1L || is.na(self@by_profile)) {
      "@by_profile must be TRUE or FALSE"
    }
  }
)

S7::method(print, mnar) <- function(x, ...) {
  cat("<mnar> verification depends on the test results and on the disease\n",
      "  odds_ratio_prior: ", prior_label(x@odds_ratio_prior),
      if (x@by_profile) ", one per profile", "\n",
      "  intercepts_prior: ", prior_label(x@intercepts_prior), "\n", sep = "")
  invisible(x)
}
