# baysesp 0.1.0

First release.

* `baysesp_data()` prepares individual or aggregated data on one or more binary
  tests, verification status and disease status; `simulate_study()` simulates
  partially verified studies with known accuracy.
* `baysesp()` fits the models in Stan, with the precompiled `{rstan}` models by
  default or with `{cmdstanr}` (`compile_models()`); `update()` refits with
  other settings.
* Three models: `saturated()` makes no assumption on the dependence between
  tests and separates identified and non-identified parameters;
  `loglinear()` and `random_effects()` assume a structure for the conditional
  dependence between tests.
* Verification at random (`mar()`) or depending on the disease (`mnar()`).
* Priors: `prior_beta()`, `elicit_beta()`, `prior_normal()`, `prior_student()`
  and `prior_dirichlet()`; profiles that are never verified require an explicit
  prior.
* Accuracy of single tests and of combination rules (`any_positive()`,
  `all_positive()`, `min_positive()`) with `test_accuracy()` and
  `rule_accuracy()`; comparisons between tests with `compare_tests()` and
  `accuracy_ratios()`.
* Quantities that the data cannot identify: `sensitivity_analysis()` over
  priors and odds ratios of verification, `accuracy_bounds()`, and
  `prior_posterior()` to compare priors and posteriors.
* `compare_models()` compares fitted models by leave-one-out cross-validation
  and posterior predictive checks.
