#' @include models.R generics.R
NULL

#' Random effects model
#'
#' Models the joint distribution of the test results with one continuous
#' latent variable \eqn{Z \sim N(0, 1)}, independent of the disease and that
#' can be read as the severity of the condition. Within each disease class the
#' tests are independent given \eqn{Z}, and they are dependent overall because
#' they share \eqn{Z}. The model needs at least two tests.
#'
#' @details
#' Let \eqn{t = (t_1, \dots, t_K)} be a response profile, with \eqn{t_k = 1}
#' for a positive result, and \eqn{d \in \{0, 1\}} the disease status. Given
#' \eqn{D = d} and \eqn{Z = z}, the test results are independent and
#' \deqn{P(T_k = 1 \mid D = d, Z = z) = \Phi(a_{dk} + b_{dk} z),}
#' where \eqn{\Phi} is the standard normal distribution function, \eqn{a_{dk}}
#' is an intercept and \eqn{b_{dk}} a loading, for each test \eqn{k} and class
#' \eqn{d}. The model has \eqn{2K} intercepts and loadings per class, plus the
#' prevalence \eqn{\pi = P(D = 1)}. This is the Gaussian random effects (probit
#' latent variable) model for conditional dependence between tests
#' (Qu et al., 1996; Dendukuri and Joseph, 2001).
#'
#' Integrating out \eqn{Z} gives the probability of each profile in class
#' \eqn{d}, \eqn{\theta_{d,t} = P(T = t \mid D = d)}:
#' \deqn{\theta_{d,t} = \int \prod_{k=1}^{K} \Phi(a_{dk} + b_{dk} z)^{t_k}\,
#' [1 - \Phi(a_{dk} + b_{dk} z)]^{1 - t_k}\,\phi(z)\,dz,}
#' where \eqn{\phi} is the standard normal density. The integral is computed
#' with 32-node Gauss-Hermite quadrature.
#'
#' Integrating out \eqn{Z} also gives the marginal probits
#' \deqn{P(T_k = 1 \mid D = d) = \Phi(m_{dk}), \qquad
#' m_{dk} = \frac{a_{dk}}{\sqrt{1 + b_{dk}^2}},}
#' so that the sensitivity of test \eqn{k} is \eqn{Se_k = \Phi(m_{1k})} and
#' its specificity is \eqn{Sp_k = 1 - \Phi(m_{0k})}. The model is sampled in
#' the marginal probits \eqn{m_{dk}}, which the data inform directly, and
#' \eqn{a_{dk} = m_{dk}\sqrt{1 + b_{dk}^2}}.
#'
#' The loadings describe the dependence between the tests within each class.
#' Given \eqn{D = d}, the covariance between tests \eqn{j} and \eqn{k} has the
#' sign of \eqn{b_{dj} b_{dk}}: a positive product means the tests tend to
#' agree beyond what the disease explains, a negative product means they tend
#' to disagree. The size of the loadings sets the strength of the dependence,
#' and a loading close to 0 makes its test almost independent of the others in
#' that class. Loadings are specific to each class, so the dependence can
#' differ between diseased and healthy subjects.
#'
#' Let \eqn{n_{t,1}} and \eqn{n_{t,0}} be the numbers of verified diseased and
#' verified healthy subjects with profile \eqn{t}, \eqn{u_t} the number of
#' unverified subjects with profile \eqn{t}, \eqn{\pi_1 = \pi} and
#' \eqn{\pi_0 = 1 - \pi}. The log-likelihood of the aggregated data is
#' \deqn{\ell = \sum_t \Big[\sum_{d} n_{t,d} \log(\pi_d\,\theta_{d,t}\,
#' \lambda_{t,d}) + u_t \log \sum_{d} \pi_d\,\theta_{d,t}\,(1 - \lambda_{t,d})\Big],}
#' where \eqn{\lambda_{t,d} = P(V = 1 \mid T = t, D = d)} is the probability of
#' verification.
#' Dependence models with partial verification of this kind are discussed by
#' Albert and Dodd (2008).
#'
#' With `verification = mar()` (missing at random) the verification probability does not depend on the
#' disease status (Rubin, 1976), so the verification terms carry no information on the
#' parameters and are dropped:
#' \deqn{\ell = \sum_t \Big[\sum_{d} n_{t,d} \log(\pi_d\,\theta_{d,t}) +
#' u_t \log \sum_{d} \pi_d\,\theta_{d,t}\Big] + \textrm{const.}}
#'
#' With `verification = mnar()` (missing not at random) the verification probability depends on the disease
#' status (Kosinski and Barnhart, 2003): \eqn{\mathrm{logit}\,\lambda_{t,d} = \alpha_t + \beta d} for every
#' profile with at least one verified subject, while \eqn{\lambda_{t,d} = 0}
#' for profiles that are never verified. The parameter \eqn{\beta} is the log
#' odds ratio of verification between diseased and healthy subjects with the
#' same profile; \eqn{\beta > 0} means that the diseased are verified more
#' often, and \eqn{\exp(-\beta)} is the odds ratio of disease between unverified
#' and verified subjects with the same profile. With `by_profile = TRUE` in
#' `verification = mnar(by_profile = TRUE)`, each profile has its own \eqn{\beta_t}.
#'
#' The data do not inform every parameter directly. The probabilities
#' \eqn{\theta_{d,t}} of the profiles that are never verified are determined
#' only by the assumed latent structure, and under MNAR the parameter
#' \eqn{\beta} is also informed only by that structure. Both depend on the
#' correctness of the model, so the results on these quantities are fragile
#' if the structure is wrong (Albert and Dodd, 2004).
#'
#' **Constraint on the signs of the loadings.** Since \eqn{Z} and \eqn{-Z}
#' have the same distribution, changing the sign of all the loadings of a class
#' gives exactly the same profile probabilities: the parameters
#' \eqn{(a_{d1}, \dots, a_{dK}, b_{d1}, \dots, b_{dK})} and
#' \eqn{(a_{d1}, \dots, a_{dK}, -b_{d1}, \dots, -b_{dK})} fit the data equally
#' well. Without a constraint the posterior would have two mirror-image modes,
#' and different chains could settle in different ones. Two constraints are
#' available:
#'
#' - `positive_loadings = TRUE` (default): all the loadings are positive. The
#'   mirror image is excluded entirely, and every pair of tests is positively
#'   dependent within each class: the tests tend to err together, as is usual
#'   for tests that measure the same condition. Negative dependence cannot be
#'   represented, and independence (all loadings 0) is on the boundary of the
#'   parameter space: when the tests are in fact conditionally independent,
#'   the posterior keeps some positive dependence. The accuracy of each test is
#'   then estimated well, but quantities driven by the dependence are not: in
#'   simulations with independent tests, the 95% interval for the
#'   sensitivity of the rule "any test positive" covered the true value in
#'   only about 75% of the datasets, because the true value lies at the edge of what the model can represent.
#'   Compare the model with `loglinear(interaction_order = 1)` in
#'   [compare_models()] before relying on such quantities.
#' - `positive_loadings = FALSE`: only the first loading of each class is
#'   positive, the others can take any sign, so any sign of dependence is
#'   allowed. This removes the mirror image only when the first loading is
#'   clearly away from 0. When it is close to 0, for example because the first
#'   test is nearly independent of the others, the mirror-image configuration
#'   (first loading near 0, the others with the opposite sign) can be reached
#'   continuously, and the posterior keeps a second mode. Chains that settle in
#'   different modes give different estimates (R-hat well above 1), and the
#'   sampler can produce divergent transitions near the boundary. On simulated
#'   data with dependent tests this happens in most fits. Use this option only
#'   when negative dependence is plausible, with several chains, and check that
#'   the chains agree.
#'
#' @param prevalence_prior Prior on the prevalence \eqn{\pi}: a single prior,
#'   either on a probability, such as `prior_beta()`, which acts on \eqn{\pi},
#'   or on the real line, such as `prior_normal()`, which acts on
#'   \eqn{\mathrm{logit}\,\pi}. Default `prior_beta(0.5, 2)`, the uniform
#'   distribution on (0, 1).
#' @param intercepts_prior Prior on the marginal probits \eqn{m_{dk}}, a
#'   real-valued prior applied to all \eqn{2K} values (\eqn{K} tests for each
#'   class) with the same prior, for example `prior_normal()`. For the diseased
#'   class it is a prior on the probit of the sensitivity, for the healthy
#'   class on the probit of one minus the specificity. Default
#'   `prior_normal(0, 1.5)`.
#' @param loadings_prior Prior on the loadings \eqn{b_{dk}}, a real-valued prior
#'   applied to all \eqn{2K} loadings with the same prior, for example
#'   `prior_normal()`. If `positive_loadings = TRUE` the prior is restricted to
#'   positive values for all loadings; if `FALSE` it is restricted to positive
#'   values only for the first loading of each class, the others are real.
#'   Default `prior_normal(0, 1)`.
#' @param positive_loadings A single logical value, `TRUE` or `FALSE`, not
#'   `NA`. If `TRUE` (default), all the loadings are positive, so the tests can
#'   only be positively dependent within each class. If `FALSE`, only the first
#'   loading of each class is positive and the others can take any sign, which
#'   can make the sampler unreliable. See Details for why a constraint is
#'   needed and how the two options behave.
#' @param verification Verification mechanism: an object of class `mar`
#'   (default, `mar()`) or `mnar()`. The priors of the verification mechanism
#'   (for example `verification = mnar(odds_ratio_prior = prior_normal(0, 1), by_profile = TRUE)`) are set in
#'   that object, not in this function.
#' @return An object of class `random_effects`.
#' @references
#' Albert, P. S. and Dodd, L. E. (2004). A cautionary note on the robustness of latent class models for estimating diagnostic error without a gold standard. *Biometrics*, 60(2), 427–435. \doi{10.1111/j.0006-341X.2004.00187.x}
#'
#' Albert, P. S. and Dodd, L. E. (2008). On estimating diagnostic accuracy from studies with multiple raters and partial gold standard evaluation. *Journal of the American Statistical Association*, 103(481), 61–73. \doi{10.1198/016214507000000329}
#'
#' Dendukuri, N. and Joseph, L. (2001). Bayesian approaches to modeling the conditional dependence between multiple diagnostic tests. *Biometrics*, 57(1), 158–167. \doi{10.1111/j.0006-341X.2001.00158.x}
#'
#' Kosinski, A. S. and Barnhart, H. X. (2003). Accounting for nonignorable verification bias in assessment of diagnostic tests. *Biometrics*, 59(1), 163–171. \doi{10.1111/1541-0420.00019}
#'
#' Qu, Y., Tan, M. and Kutner, M. H. (1996). Random effects models in latent class analysis for evaluating accuracy of diagnostic tests. *Biometrics*, 52(3), 797–810. \doi{10.2307/2533043}
#'
#' Rubin, D. B. (1976). Inference and missing data. *Biometrika*, 63(3), 581–592. \doi{10.1093/biomet/63.3.581}
#' @examples
#' random_effects()
#' random_effects(positive_loadings = FALSE)
#' @export
random_effects <- S7::new_class(
  "random_effects",
  parent = accuracy_model,
  properties = list(
    prevalence_prior = prior_distribution,
    intercepts_prior = prior_distribution,
    loadings_prior = prior_distribution,
    positive_loadings = S7::class_logical
  ),
  constructor = function(prevalence_prior = prior_beta(0.5, 2),
                         intercepts_prior = prior_normal(0, 1.5),
                         loadings_prior = prior_normal(0, 1),
                         positive_loadings = TRUE,
                         verification = mar()) {
    S7::new_object(
      S7::S7_object(),
      verification = verification,
      prevalence_prior = prevalence_prior,
      intercepts_prior = intercepts_prior,
      loadings_prior = loadings_prior,
      positive_loadings = positive_loadings
    )
  },
  validator = function(self) {
    targets <- c(prevalence_prior = "probability", intercepts_prior = "real",
                 loadings_prior = "real")
    for (name in names(targets)) {
      problem <- tryCatch(check_prior(S7::prop(self, name), targets[[name]], name),
                          error = conditionMessage)
      if (is.character(problem)) return(problem)
    }
    if (length(self@positive_loadings) != 1L || is.na(self@positive_loadings)) {
      "@positive_loadings must be TRUE or FALSE"
    }
  }
)

S7::method(print_generic, random_effects) <- function(x, ...) {
  cat("<random_effects> random effects model, Gauss-Hermite quadrature\n",
      "  prevalence_prior: ", prior_label(x@prevalence_prior), "\n",
      "  intercepts_prior: ", prior_label(x@intercepts_prior), "\n",
      "  loadings_prior:   ", prior_label(x@loadings_prior),
      if (x@positive_loadings) ", all loadings positive" else ", first loading positive",
      "\n",
      "  verification:     ", sub("^baysesp::", "", class(x@verification)[1]),
      "\n", sep = "")
  invisible(x)
}
