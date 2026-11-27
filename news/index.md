# Changelog

## baysesp 0.1.0

First release.

- [`baysesp_data()`](https://giovannitinervia9.github.io/baysesp/reference/baysesp_data.md)
  prepares individual or aggregated data on one or more binary tests,
  verification status and disease status;
  [`simulate_study()`](https://giovannitinervia9.github.io/baysesp/reference/simulate_study.md)
  simulates partially verified studies with known accuracy.
- [`baysesp()`](https://giovannitinervia9.github.io/baysesp/reference/baysesp.md)
  fits the models in Stan, with the precompiled
  [rstan](https://mc-stan.org/rstan/) models by default or with
  `{cmdstanr}`
  ([`compile_models()`](https://giovannitinervia9.github.io/baysesp/reference/compile_models.md));
  [`update()`](https://rdrr.io/r/stats/update.html) refits with other
  settings.
- Three models:
  [`saturated()`](https://giovannitinervia9.github.io/baysesp/reference/saturated.md)
  makes no assumption on the dependence between tests and separates
  identified and non-identified parameters;
  [`loglinear()`](https://giovannitinervia9.github.io/baysesp/reference/loglinear.md)
  and
  [`random_effects()`](https://giovannitinervia9.github.io/baysesp/reference/random_effects.md)
  assume a structure for the conditional dependence between tests.
- Verification at random
  ([`mar()`](https://giovannitinervia9.github.io/baysesp/reference/mar.md))
  or depending on the disease
  ([`mnar()`](https://giovannitinervia9.github.io/baysesp/reference/mnar.md)).
- Priors:
  [`prior_beta()`](https://giovannitinervia9.github.io/baysesp/reference/prior_beta.md),
  [`elicit_beta()`](https://giovannitinervia9.github.io/baysesp/reference/elicit_beta.md),
  [`prior_normal()`](https://giovannitinervia9.github.io/baysesp/reference/prior_normal.md),
  [`prior_student()`](https://giovannitinervia9.github.io/baysesp/reference/prior_student.md)
  and
  [`prior_dirichlet()`](https://giovannitinervia9.github.io/baysesp/reference/prior_dirichlet.md);
  profiles that are never verified require an explicit prior.
- Accuracy of single tests and of combination rules
  ([`any_positive()`](https://giovannitinervia9.github.io/baysesp/reference/any_positive.md),
  [`all_positive()`](https://giovannitinervia9.github.io/baysesp/reference/all_positive.md),
  [`min_positive()`](https://giovannitinervia9.github.io/baysesp/reference/min_positive.md))
  with
  [`test_accuracy()`](https://giovannitinervia9.github.io/baysesp/reference/test_accuracy.md)
  and
  [`rule_accuracy()`](https://giovannitinervia9.github.io/baysesp/reference/rule_accuracy.md);
  comparisons between tests with
  [`compare_tests()`](https://giovannitinervia9.github.io/baysesp/reference/compare_tests.md)
  and
  [`accuracy_ratios()`](https://giovannitinervia9.github.io/baysesp/reference/accuracy_ratios.md).
- Quantities that the data cannot identify:
  [`sensitivity_analysis()`](https://giovannitinervia9.github.io/baysesp/reference/sensitivity_analysis.md)
  over priors and odds ratios of verification,
  [`accuracy_bounds()`](https://giovannitinervia9.github.io/baysesp/reference/accuracy_bounds.md),
  and
  [`prior_posterior()`](https://giovannitinervia9.github.io/baysesp/reference/prior_posterior.md)
  to compare priors and posteriors.
- [`compare_models()`](https://giovannitinervia9.github.io/baysesp/reference/compare_models.md)
  compares fitted models by leave-one-out cross-validation and posterior
  predictive checks.
