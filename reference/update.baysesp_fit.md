# Fit a different model to the same data

Fits a new model to the data of a model already fitted with
[`baysesp()`](https://giovannitinervia9.github.io/baysesp/reference/baysesp.md),
with the same backend and the same sampling settings (number of chains,
warmup and sampling iterations, seed, and the further arguments given to
[`baysesp()`](https://giovannitinervia9.github.io/baysesp/reference/baysesp.md)).
Only the model changes: its class, its priors or its verification
mechanism. This is the usual way to fit the alternative models given to
[`compare_models()`](https://giovannitinervia9.github.io/baysesp/reference/compare_models.md).

## Usage

``` r
# S3 method for class 'baysesp_fit'
update(object, model = object@model, ...)
```

## Arguments

- object:

  A model fitted with
  [`baysesp()`](https://giovannitinervia9.github.io/baysesp/reference/baysesp.md).

- model:

  The new model: an object created by
  [`saturated()`](https://giovannitinervia9.github.io/baysesp/reference/saturated.md),
  [`loglinear()`](https://giovannitinervia9.github.io/baysesp/reference/loglinear.md)
  or
  [`random_effects()`](https://giovannitinervia9.github.io/baysesp/reference/random_effects.md),
  with its priors and its verification mechanism (`verification = mar()`
  or `verification = mnar()`). Default `object@model`, the model of
  `object`, which fits the same model again.

- ...:

  Sampling settings that replace those of `object`, with the names of
  the arguments of
  [`baysesp()`](https://giovannitinervia9.github.io/baysesp/reference/baysesp.md):
  `backend`, `chains`, `iter_warmup`, `iter_sampling`, `seed`, or
  arguments for the backend such as `cores`. Settings that are not given
  are those used for `object`; the seed, if `object` was fitted with
  `seed = NULL`, is the one drawn at that time.

## Value

A new `baysesp_fit` object, as returned by
[`baysesp()`](https://giovannitinervia9.github.io/baysesp/reference/baysesp.md).

## See also

[`baysesp()`](https://giovannitinervia9.github.io/baysesp/reference/baysesp.md),
[`compare_models()`](https://giovannitinervia9.github.io/baysesp/reference/compare_models.md)

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
fit_ci <- update(fit, model = loglinear(interaction_order = 1))
fit_ci
#> <baysesp_fit> loglinear model, MAR verification, 4 chains x 1000 draws (rstan)
#> Convergence: max R-hat 1.002, min bulk ESS 4426, min tail ESS 2489, 0 divergent transitions
#> 
#>    quantity  mean median lower upper identified
#>  prevalence 0.205  0.205 0.180 0.232       TRUE
#>       Se[A] 0.885  0.886 0.831 0.929       TRUE
#>       Sp[A] 0.849  0.849 0.823 0.874       TRUE
#>      PPV[A] 0.602  0.602 0.546 0.657       TRUE
#>      NPV[A] 0.966  0.967 0.950 0.980       TRUE
#>       Se[B] 0.788  0.789 0.726 0.844       TRUE
#>       Sp[B] 0.946  0.947 0.929 0.962       TRUE
#>      PPV[B] 0.791  0.792 0.731 0.847       TRUE
#>      NPV[B] 0.945  0.946 0.927 0.962       TRUE
#> 
#> Odds ratio of verification fixed at 1 (MAR).
# }
```
