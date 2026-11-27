#' @include priors.R generics.R
NULL

#' Verification mechanisms
#'
#' `verification_mechanism` is the abstract parent class of the assumptions on
#' how subjects are selected for verification by the reference standard. It
#' cannot be instantiated: the concrete assumptions are the missing at random
#' assumption, in which verification depends only on the test results, and the
#' missing not at random assumption, in which verification also depends on the
#' disease status (Rubin, 1976; Little and Rubin, 2002). Each concrete assumption has its own page.
#'
#' @references
#' Little, R. J. A. and Rubin, D. B. (2002). *Statistical Analysis with Missing Data* (2nd ed.). Wiley.
#'
#' Rubin, D. B. (1976). Inference and missing data. *Biometrika*, 63(3), 581–592. \doi{10.1093/biomet/63.3.581}
#' @seealso [mar()], [mnar()]
#' @examples
#' S7::S7_inherits(mar(), verification_mechanism)
#' @export
verification_mechanism <- S7::new_class("verification_mechanism", abstract = TRUE)

#' Missing at random verification
#'
#' `mar()` is the assumption that verification depends only on the test results.
#' Once the response profile is known, the disease status does not change the
#' probability of being verified. The function has no arguments.
#'
#' @details
#' Let \eqn{V = 1} indicate that a subject is verified, \eqn{T} be the response
#' profile and \eqn{D} the disease status, and let
#' \eqn{\lambda_{t,d} = P(V = 1 \mid T = t, D = d)}. The assumption is
#' \deqn{\lambda_{t,1} = \lambda_{t,0} = \lambda_t \quad \textrm{for every profile } t,}
#' so the probability of verification depends on the profile only.
#'
#' Equivalently, among subjects with the same profile \eqn{t}, verified and
#' unverified subjects have the same proportion of diseased:
#' \deqn{P(D = 1 \mid T = t, V = 1) = P(D = 1 \mid T = t, V = 0) = P(D = 1 \mid T = t).}
#' Verification therefore carries no information on the disease beyond the
#' test results.
#'
#' Under this assumption the verification process is ignorable and can be left
#' out of the likelihood (Rubin, 1976; Little and Rubin, 2002). Let \eqn{n_{t,d}} be the number of verified subjects
#' with profile \eqn{t} and disease status \eqn{d}, \eqn{u_t} the number of
#' unverified subjects with profile \eqn{t}, \eqn{\pi = P(D = 1)} the
#' prevalence, \eqn{\pi_1 = \pi}, \eqn{\pi_0 = 1 - \pi}, and
#' \eqn{\theta_{d,t} = P(T = t \mid D = d)}. The log-likelihood of the observed
#' counts is
#' \deqn{\ell = \sum_t \Big[\sum_d n_{t,d} \log(\pi_d\,\theta_{d,t}\,\lambda_{t,d}) + u_t \log \sum_d \pi_d\,\theta_{d,t}\,(1 - \lambda_{t,d})\Big].}
#' Setting \eqn{\lambda_{t,d} = \lambda_t} splits it into two parts,
#' \deqn{\ell = \ell_D + \ell_V,}
#' \deqn{\ell_D = \sum_t \Big[\sum_d n_{t,d} \log(\pi_d\,\theta_{d,t}) + u_t \log \sum_d \pi_d\,\theta_{d,t}\Big], \qquad \ell_V = \sum_t \Big[(n_{t,1} + n_{t,0}) \log \lambda_t + u_t \log(1 - \lambda_t)\Big].}
#' The verification probabilities \eqn{\lambda_t} appear only in \eqn{\ell_V}.
#' Inference on the prevalence and on the test distributions is based on
#' \eqn{\ell_D}, and the unverified subjects contribute to it through the
#' totals \eqn{u_t} (Begg and Greenes, 1983; Zhou, 1993).
#'
#' MAR is plausible when the choice of who undergoes the reference standard is
#' made from the test results alone, for example when a protocol prescribes the
#' reference test for every subject with a given test pattern, and no other
#' information related to the disease affects the choice. It is doubtful when
#' the decision to verify uses information on the disease that the tests do not
#' capture, such as symptoms, medical history or clinical suspicion, or when,
#' within the same profile, subjects who are more likely to be diseased are
#' verified more often. In those cases the proportion of diseased among verified
#' subjects differs from that among unverified subjects with the same profile,
#' and the estimates obtained under MAR are biased. The data do not contain
#' information on whether MAR holds, so its plausibility must be judged on
#' substantive grounds. The alternative in which verification also depends on
#' the disease status is missing not at random.
#'
#' Profiles that are never verified are compatible with MAR, but for them the
#' data give no information on the proportion of diseased: there are no
#' verified subjects from which to estimate it. What is concluded about it
#' comes entirely from an assumption: a prior on that proportion (Buzoianu and Kadane, 2008; Gustafson, 2015), or a
#' structure imposed on the dependence between the tests.
#'
#' @return An object of class `mar`.
#' @references
#' Begg, C. B. and Greenes, R. A. (1983). Assessment of diagnostic tests when disease verification is subject to selection bias. *Biometrics*, 39(1), 207–215. \doi{10.2307/2530820}
#'
#' Buzoianu, M. and Kadane, J. B. (2008). Adjusting for verification bias in diagnostic test evaluation: A Bayesian approach. *Statistics in Medicine*, 27(13), 2453–2473. \doi{10.1002/sim.3099}
#'
#' Gustafson, P. (2015). *Bayesian Inference for Partially Identified Models: Exploring the Limits of Limited Data*. Chapman & Hall/CRC, Boca Raton, FL.
#'
#' Little, R. J. A. and Rubin, D. B. (2002). *Statistical Analysis with Missing Data* (2nd ed.). Wiley.
#'
#' Rubin, D. B. (1976). Inference and missing data. *Biometrika*, 63(3), 581–592. \doi{10.1093/biomet/63.3.581}
#'
#' Zhou, X.-H. (1993). Maximum likelihood estimators of sensitivity and specificity corrected for verification bias. *Communications in Statistics -- Theory and Methods*, 22(11), 3177–3198. \doi{10.1080/03610929308831209}
#' @seealso [mnar()], [verification_mechanism]
#' @examples
#' mar()
#' @export
mar <- S7::new_class("mar", parent = verification_mechanism)

S7::method(print_generic, mar) <- function(x, ...) {
  cat("<mar> verification depends only on the test results\n")
  invisible(x)
}

#' Missing not at random verification
#'
#' `mnar()` is the assumption that verification depends on the test results and
#' also on the disease status (Kosinski and Barnhart, 2003a): among subjects with the same profile, diseased
#' and healthy subjects can be verified with different probabilities. Its
#' arguments are the prior on the log odds ratio of verification, the prior on
#' the verification intercepts and the choice of whether the odds ratio varies
#' across profiles.
#'
#' @details
#' Let \eqn{V = 1} indicate that a subject is verified, \eqn{T} be the response
#' profile, \eqn{D} the disease status, and
#' \eqn{\lambda_{t,d} = P(V = 1 \mid T = t, D = d)}. The model is
#' \deqn{\mathrm{logit}\,\lambda_{t,d} = \alpha_t + \beta d, \qquad d \in \{0, 1\}.}
#' For each profile \eqn{t} with at least one verified subject, \eqn{\alpha_t}
#' is the log odds of verification for healthy subjects with profile \eqn{t},
#' and \eqn{\alpha_t + \beta} is the log odds of verification for diseased
#' subjects with the same profile. Profiles with no verified subject have
#' \eqn{\lambda_{t,d} = 0} and no intercept; in the log-linear and random
#' effects models `baysesp()` warns when such a profile has subjects, because
#' for it verification is then not missing not at random.
#'
#' \eqn{\beta} is the log odds ratio of verification between diseased and
#' healthy subjects with the same profile:
#' \deqn{\beta = \log \frac{\lambda_{t,1}\,/\,(1 - \lambda_{t,1})}{\lambda_{t,0}\,/\,(1 - \lambda_{t,0})}.}
#' Equivalently, \eqn{\exp(-\beta)} is the odds ratio of disease between
#' unverified and verified subjects with the same profile. When \eqn{\beta = 0}
#' verification does not depend on the disease status given the profile, which
#' is the missing at random assumption. A positive \eqn{\beta}
#' means that diseased subjects are verified more often than healthy subjects
#' with the same profile, so that unverified subjects are less often diseased
#' than verified ones with the same profile; a negative \eqn{\beta} has the
#' opposite meaning.
#'
#' With \eqn{\pi = P(D = 1)} and \eqn{\theta_{d,t} = P(T = t \mid D = d)} the
#' disease parameters, and \eqn{n_{t,d}}, \eqn{u_t}, \eqn{\pi_1 = \pi},
#' \eqn{\pi_0 = 1 - \pi} as in the counts of the observed data, the
#' log-likelihood is
#' \deqn{\ell = \sum_t \Big[\sum_d n_{t,d} \log(\pi_d\,\theta_{d,t}\,\lambda_{t,d}) + u_t \log \sum_d \pi_d\,\theta_{d,t}\,(1 - \lambda_{t,d})\Big].}
#'
#' With a saturated disease model, which places no restriction on
#' \eqn{\theta_d}, the data identify the proportion of diseased among verified
#' subjects of each profile, but not \eqn{\beta}. The proportion of diseased
#' among unverified subjects with a verified profile \eqn{t} is
#' \eqn{\mathrm{expit}(\mathrm{logit}\,\phi_t - \beta)}, where \eqn{\phi_t} is
#' the proportion among verified subjects, so \eqn{\beta} links the two groups
#' but the data say nothing about it: its posterior equals its prior (Gustafson, 2015). The prior
#' on \eqn{\beta} is then an explicit assumption, and varying it, or fixing
#' \eqn{\beta} at a range of values, shows how much the results depend on it (Buzoianu and Kadane, 2008; Kosinski and Barnhart, 2003b).
#'
#' With the log-linear (Baker, 1995) and random-effects (Dendukuri and Joseph, 2001; Qu et al., 1996) disease models, \eqn{\beta} is
#' identified through the restrictions that these models place on the disease
#' distribution. This information comes from the assumed structure, not from
#' the data alone, and it is fragile if the structure is wrong (Albert and Dodd, 2004).
#'
#' @param odds_ratio_prior Prior on \eqn{\beta}, the log odds ratio of
#'   verification between diseased and healthy subjects with the same profile.
#'   It must act on the real line, such as `prior_normal()` or
#'   `prior_student()`; a prior on a probability, such as `prior_beta()`, is
#'   not accepted. The default `prior_normal(0, 1)` centers \eqn{\beta} at 0,
#'   which is missing at random, and gives 95% of its prior mass to odds
#'   ratios between about 0.14 and 7.1.
#' @param intercepts_prior Prior on the verification intercepts \eqn{\alpha_t},
#'   one for each profile with at least one verified subject. It must act on
#'   the real line, as for `odds_ratio_prior`. It is used only by the log-linear
#'   and the random-effects models; the saturated model does not use it. The
#'   default is `prior_normal(0, 2.5)`. In profiles where almost all subjects
#'   are diseased, the data say little about the verification probability of
#'   the healthy, \eqn{\mathrm{expit}(\alpha_t)}, and the prior pulls it
#'   towards 0.5; the shared \eqn{\beta} then absorbs the difference. As a
#'   result the estimated odds ratio depends on this prior: in simulations
#'   with verification probabilities of 0.8, `prior_normal(0, 1.5)` gave a
#'   posterior median of the odds ratio of about 2.4 when the true value was
#'   1, the default about 1.7. The sensitivities and specificities changed by
#'   less than 0.005. Read the odds ratio of the structured models together
#'   with [sensitivity_analysis()].
#' @param by_profile A single logical value, `TRUE` or `FALSE`, with default
#'   `FALSE`. If `FALSE`, one \eqn{\beta} is shared by all profiles with at
#'   least one verified subject, so the odds ratio of verification is the same
#'   for every profile. If `TRUE`, each of these profiles has its own
#'   \eqn{\beta_t}, all with prior `odds_ratio_prior`, so the odds ratio may
#'   differ across profiles. Profiles with no verified subject have no
#'   \eqn{\beta}.
#' @return An object of class `mnar`.
#' @references
#' Albert, P. S. and Dodd, L. E. (2004). A cautionary note on the robustness of latent class models for estimating diagnostic error without a gold standard. *Biometrics*, 60(2), 427–435. \doi{10.1111/j.0006-341X.2004.00187.x}
#'
#' Baker, S. G. (1995). Evaluating multiple diagnostic tests with partial verification. *Biometrics*, 51(1), 330–337. \doi{10.2307/2533339}
#'
#' Buzoianu, M. and Kadane, J. B. (2008). Adjusting for verification bias in diagnostic test evaluation: A Bayesian approach. *Statistics in Medicine*, 27(13), 2453–2473. \doi{10.1002/sim.3099}
#'
#' Dendukuri, N. and Joseph, L. (2001). Bayesian approaches to modeling the conditional dependence between multiple diagnostic tests. *Biometrics*, 57(1), 158–167. \doi{10.1111/j.0006-341X.2001.00158.x}
#'
#' Gustafson, P. (2015). *Bayesian Inference for Partially Identified Models: Exploring the Limits of Limited Data*. Chapman & Hall/CRC, Boca Raton, FL.
#'
#' Kosinski, A. S. and Barnhart, H. X. (2003a). Accounting for nonignorable verification bias in assessment of diagnostic tests. *Biometrics*, 59(1), 163–171. \doi{10.1111/1541-0420.00019}
#'
#' Kosinski, A. S. and Barnhart, H. X. (2003b). A global sensitivity analysis of performance of a medical diagnostic test when verification bias is present. *Statistics in Medicine*, 22(17), 2711–2721. \doi{10.1002/sim.1517}
#'
#' Qu, Y., Tan, M. and Kutner, M. H. (1996). Random effects models in latent class analysis for evaluating accuracy of diagnostic tests. *Biometrics*, 52(3), 797–810. \doi{10.2307/2533043}
#' @seealso [mar()], [verification_mechanism], [sensitivity_analysis()]
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

S7::method(print_generic, mnar) <- function(x, ...) {
  cat("<mnar> verification depends on the test results and on the disease\n",
      "  odds_ratio_prior: ", prior_label(x@odds_ratio_prior),
      if (x@by_profile) ", one per profile", "\n",
      "  intercepts_prior: ", prior_label(x@intercepts_prior), "\n", sep = "")
  invisible(x)
}
