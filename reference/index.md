# Package index

## Data

- [`baysesp_data()`](https://giovannitinervia9.github.io/baysesp/reference/baysesp_data.md)
  : Study data
- [`simulate_study()`](https://giovannitinervia9.github.io/baysesp/reference/simulate_study.md)
  : Simulate a study with partial verification
- [`independence()`](https://giovannitinervia9.github.io/baysesp/reference/independence.md)
  : Independent tests in simulated data
- [`latent_dependence()`](https://giovannitinervia9.github.io/baysesp/reference/latent_dependence.md)
  : Tests dependent through a shared latent variable in simulated data
- [`probit_dependence()`](https://giovannitinervia9.github.io/baysesp/reference/probit_dependence.md)
  : Tests with correlated latent variables in simulated data

## Models and fitting

- [`baysesp()`](https://giovannitinervia9.github.io/baysesp/reference/baysesp.md)
  : Bayesian estimation of the accuracy of diagnostic tests
- [`saturated()`](https://giovannitinervia9.github.io/baysesp/reference/saturated.md)
  : Saturated model
- [`loglinear()`](https://giovannitinervia9.github.io/baysesp/reference/loglinear.md)
  : Log-linear model of the test profiles
- [`random_effects()`](https://giovannitinervia9.github.io/baysesp/reference/random_effects.md)
  : Random effects model
- [`accuracy_model()`](https://giovannitinervia9.github.io/baysesp/reference/accuracy_model.md)
  : Accuracy models
- [`compile_models()`](https://giovannitinervia9.github.io/baysesp/reference/compile_models.md)
  : Compile the Stan programs for the cmdstanr backend
- [`summary(`*`<baysesp_fit>`*`)`](https://giovannitinervia9.github.io/baysesp/reference/summary.baysesp_fit.md)
  : Summary of a fitted model
- [`update(`*`<baysesp_fit>`*`)`](https://giovannitinervia9.github.io/baysesp/reference/update.baysesp_fit.md)
  : Fit a different model to the same data

## Verification

- [`mar()`](https://giovannitinervia9.github.io/baysesp/reference/mar.md)
  : Missing at random verification
- [`mnar()`](https://giovannitinervia9.github.io/baysesp/reference/mnar.md)
  : Missing not at random verification
- [`verification_mechanism()`](https://giovannitinervia9.github.io/baysesp/reference/verification_mechanism.md)
  : Verification mechanisms

## Priors

- [`prior_beta()`](https://giovannitinervia9.github.io/baysesp/reference/prior_beta.md)
  : Beta prior with mean and precision
- [`elicit_beta()`](https://giovannitinervia9.github.io/baysesp/reference/elicit_beta.md)
  : Beta prior elicited from a mean and a plausible upper value
- [`prior_normal()`](https://giovannitinervia9.github.io/baysesp/reference/prior_normal.md)
  : Normal prior
- [`prior_student()`](https://giovannitinervia9.github.io/baysesp/reference/prior_student.md)
  : Student t prior
- [`prior_dirichlet()`](https://giovannitinervia9.github.io/baysesp/reference/prior_dirichlet.md)
  : Dirichlet prior
- [`prior_distribution()`](https://giovannitinervia9.github.io/baysesp/reference/prior_distribution.md)
  : Prior distributions

## Diagnostic rules

- [`single_test()`](https://giovannitinervia9.github.io/baysesp/reference/single_test.md)
  : Single test rule
- [`any_positive()`](https://giovannitinervia9.github.io/baysesp/reference/any_positive.md)
  : Any positive rule (OR rule)
- [`all_positive()`](https://giovannitinervia9.github.io/baysesp/reference/all_positive.md)
  : All positive rule (AND rule)
- [`min_positive()`](https://giovannitinervia9.github.io/baysesp/reference/min_positive.md)
  : Minimum positive rule
- [`test_rule()`](https://giovannitinervia9.github.io/baysesp/reference/test_rule.md)
  : Diagnostic rule class

## Accuracy of tests and rules

- [`test_accuracy()`](https://giovannitinervia9.github.io/baysesp/reference/test_accuracy.md)
  : Accuracy of the individual tests
- [`rule_accuracy()`](https://giovannitinervia9.github.io/baysesp/reference/rule_accuracy.md)
  : Accuracy of diagnostic rules
- [`compare_tests()`](https://giovannitinervia9.github.io/baysesp/reference/compare_tests.md)
  : Compare two tests or rules

## Never-verified profiles and sensitivity analysis

- [`accuracy_bounds()`](https://giovannitinervia9.github.io/baysesp/reference/accuracy_bounds.md)
  : Bounds for the absolute accuracy when some profiles are never
  verified
- [`accuracy_ratios()`](https://giovannitinervia9.github.io/baysesp/reference/accuracy_ratios.md)
  : Ratios between the accuracy of tests when some profiles are never
  verified
- [`sensitivity_analysis()`](https://giovannitinervia9.github.io/baysesp/reference/sensitivity_analysis.md)
  : Sensitivity analysis for the assumptions the data cannot check
- [`prior_posterior()`](https://giovannitinervia9.github.io/baysesp/reference/prior_posterior.md)
  : Compare prior and posterior
- [`plot(`*`<prior_posterior>`*`)`](https://giovannitinervia9.github.io/baysesp/reference/plot.prior_posterior.md)
  : Plot the comparison of prior and posterior

## Model comparison

- [`compare_models()`](https://giovannitinervia9.github.io/baysesp/reference/compare_models.md)
  : Compare models fitted to the same data
- [`plot(`*`<model_comparison>`*`)`](https://giovannitinervia9.github.io/baysesp/reference/plot.model_comparison.md)
  : Plot the posterior predictive checks of a model comparison
