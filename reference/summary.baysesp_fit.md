# Summary of a fitted model

Collects in one object what describes a model fitted with
[`baysesp()`](https://giovannitinervia9.github.io/baysesp/reference/baysesp.md):
the model and its verification mechanism, the priors, the profiles of
test results with no verified subject, the convergence diagnostics of
the sampler and the posterior summaries of the accuracy of the tests.

## Usage

``` r
# S3 method for class 'baysesp_fit'
summary(object, level = 0.95, ...)
```

## Arguments

- object:

  A model fitted with
  [`baysesp()`](https://giovannitinervia9.github.io/baysesp/reference/baysesp.md).

- level:

  Probability of the credible intervals in `accuracy`, a single number
  in (0, 1). Default `0.95`.

- ...:

  Not used.

## Value

A list of class `fit_summary`, with a
[`print()`](https://rdrr.io/r/base/print.html) method, with the
elements:

- `model`: a description of the model, such as
  `"log-linear model of order 1"`.

- `verification`: `"MAR"` or `"MNAR"`, followed by
  `"(one odds ratio per profile)"` with `by_profile = TRUE`.

- `priors`: a data frame with the name of each argument that holds a
  prior (`argument`) and the prior (`prior`). The prior on the
  proportion of diseased in the never-verified profiles has one row per
  never-verified profile, named after it, such as
  `"unverified_prevalence_prior (0,1)"`. The priors of the verification
  mechanism are included when they are used, with the name of the
  argument of
  [`mnar()`](https://giovannitinervia9.github.io/baysesp/reference/mnar.md):
  `"mnar(odds_ratio_prior)"` with `verification = mnar()`, and
  `"mnar(intercepts_prior)"` with `verification = mnar()` in the
  log-linear and random effects models.

- `never_verified`: the never-verified profiles, written as `"(0,1)"`,
  or a character vector of length 0.

- `diagnostics`: a one-row data frame with `chains`, `draws` (per
  chain), `backend`, `max_rhat`, `min_ess_bulk`, `min_ess_tail` and
  `divergent` (number of divergent transitions; `NA` if not available).

- `accuracy`: the result of
  [`test_accuracy()`](https://giovannitinervia9.github.io/baysesp/reference/test_accuracy.md),
  with the odds ratio of verification added under
  `verification = mnar()`.

- `level`: the value of `level`.

## Details

The convergence diagnostics are computed over all the variables saved in
the draws (variables that are constant, such as the verified fraction of
a never-verified profile, are skipped): the largest split-\\\hat{R}\\,
which should be below 1.01, the smallest bulk and tail effective sample
sizes, which should be at least 100 per chain (Vehtari et al., 2021),
and the number of divergent transitions after warmup, which should be 0:
divergences mean that the sampler could not explore some region of the
posterior, and the draws may be biased.

## References

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Bürkner, P.-C.
(2021). Rank-normalization, folding, and localization: An improved
\\\hat{R}\\ for assessing convergence of MCMC (with discussion).
*Bayesian Analysis*, 16(2), 667–718.
[doi:10.1214/20-BA1221](https://doi.org/10.1214/20-BA1221)

## See also

[`baysesp()`](https://giovannitinervia9.github.io/baysesp/reference/baysesp.md),
[`test_accuracy()`](https://giovannitinervia9.github.io/baysesp/reference/test_accuracy.md)

## Examples

``` r
# \donttest{
set.seed(1)
d <- simulate_study(n = 1000, prevalence = 0.2,
                    sensitivity = c(A = 0.9, B = 0.8),
                    specificity = c(A = 0.85, B = 0.95),
                    verification_prob = c(0.2, 0.9, 0.9, 1))
d <- baysesp_data(d, tests = c("A", "B"), verified = "verified", disease = "disease")
fit <- baysesp(d, seed = 1)
summary(fit)
#> Model:          saturated model, transparent parametrization
#> Verification:   MAR
#> Never verified: none
#> 
#> Priors:
#>   frequency_prior           Dirichlet(concentration = 1)
#>   verified_prior            Beta(mean = 0.5, precision = 2)
#>   profile_prevalence_prior  Beta(mean = 0.5, precision = 2)
#> 
#> Sampling: 4 chains x 1000 draws (rstan)
#>   max R-hat 1.002, min bulk ESS 4568, min tail ESS 1724, divergent transitions 0
#> 
#> Accuracy (95% credible intervals):
#>    quantity  mean median lower upper identified
#>  prevalence 0.205  0.205 0.180 0.232       TRUE
#>       Se[A] 0.885  0.887 0.823 0.932       TRUE
#>       Sp[A] 0.849  0.849 0.823 0.873       TRUE
#>      PPV[A] 0.602  0.602 0.545 0.659       TRUE
#>      NPV[A] 0.966  0.967 0.946 0.981       TRUE
#>       Se[B] 0.788  0.789 0.719 0.851       TRUE
#>       Sp[B] 0.947  0.947 0.930 0.961       TRUE
#>      PPV[B] 0.792  0.793 0.734 0.845       TRUE
#>      NPV[B] 0.945  0.946 0.924 0.963       TRUE
# }
```
