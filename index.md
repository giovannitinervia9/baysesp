# baysesp

[baysesp](https://github.com/giovannitinervia9/baysesp) provides
Bayesian estimation of the accuracy (sensitivity, specificity and
predictive values) of one or more binary diagnostic tests performed on
the same subjects, when the disease status is ascertained by the
reference standard only for a subset of them (*partial verification*).

Estimating accuracy from the verified subjects alone is biased whenever
the decision to verify depends on the test results, as it usually does
in screening studies.
[baysesp](https://github.com/giovannitinervia9/baysesp) corrects this
*verification bias* with a joint model for the tests, the disease and
the verification process, written in [Stan](https://mc-stan.org).

## Features

- Joint modeling of **K binary tests**. The saturated model makes no
  assumption on how the tests depend on each other; log-linear and
  random effects models assume a structure for that dependence.
- Accuracy of **combination rules**: any test positive, all positive, at
  least *k* positive, or a custom rule.
- **Comparisons between tests**: differences and ratios of sensitivities
  and specificities, with their posterior distributions.
- Verification that depends only on the test results (missing at random)
  or also on the disease status (missing not at random).
- **Profiles that are never verified**, for example when only subjects
  with a positive test receive the reference standard: the missing
  information is an explicit prior, with bounds and ratios that hold
  whatever that prior is.
- **Sensitivity analyses** over the assumptions that the data cannot
  check, without fitting the model again, and prior-posterior
  comparisons.
- Two Stan backends: [rstan](https://mc-stan.org/rstan/) (default,
  models precompiled at installation) and `{cmdstanr}` (optional).

## Installation

You can install the development version of
[baysesp](https://github.com/giovannitinervia9/baysesp) from
[GitHub](https://github.com/giovannitinervia9/baysesp) with:

``` r

# install.packages("remotes")
remotes::install_github("giovannitinervia9/baysesp")
```

Installing the package compiles its Stan models and takes a few minutes.

To use the optional `{cmdstanr}` backend, install [`{cmdstanr}` and
CmdStan](https://mc-stan.org/cmdstanr/) and compile the models once with
[`compile_models()`](https://giovannitinervia9.github.io/baysesp/reference/compile_models.md):

``` r

baysesp::compile_models()
```

## Example

Two tests, `A` and `B`, are performed on 3000 simulated subjects; those
negative to both tests are rarely verified.

``` r

library(baysesp)

set.seed(2026)
study <- simulate_study(
  n = 3000,
  prevalence = 0.15,
  sensitivity = c(A = 0.85, B = 0.70),
  specificity = c(A = 0.90, B = 0.95),
  verification_prob = c(0.1, 0.8, 0.8, 1)
)

# count the subjects of each combination of test results
d <- baysesp_data(
  data = study,
  tests = c("A", "B"),
  verified = "verified",
  disease = "disease"
)

# saturated model, verification that depends only on the test results
fit <- baysesp(data = d, model = saturated(verification = mar()), seed = 1)

test_accuracy(fit = fit)                                  # Se, Sp, PPV, NPV
rule_accuracy(fit = fit, rules = list(any_positive(), all_positive()))
compare_tests(fit = fit, first = "A", second = "B")
plot(fit)
```

The [package website](https://giovannitinervia9.github.io/baysesp/) has
the documentation of every function and a vignette that goes through all
the features of the package,
[`vignette("baysesp")`](https://giovannitinervia9.github.io/baysesp/articles/baysesp.md).

## Status

The package is under active development and its interface may still
change.
