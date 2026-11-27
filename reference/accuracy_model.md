# Accuracy models

`accuracy_model` is the abstract parent class of the models for the
joint distribution of the test results and of the disease status. It
cannot be created directly: the concrete models are
[`saturated()`](https://giovannitinervia9.github.io/baysesp/reference/saturated.md),
[`loglinear()`](https://giovannitinervia9.github.io/baysesp/reference/loglinear.md)
and
[`random_effects()`](https://giovannitinervia9.github.io/baysesp/reference/random_effects.md),
each with its own constructor.

## Usage

``` r
accuracy_model(verification = verification_mechanism())
```

## Arguments

- verification:

  Verification mechanism: an object created with
  [`mar()`](https://giovannitinervia9.github.io/baysesp/reference/mar.md)
  or
  [`mnar()`](https://giovannitinervia9.github.io/baysesp/reference/mnar.md).

## Details

Every model has the property `verification`, an object of class
`verification_mechanism`
([`mar()`](https://giovannitinervia9.github.io/baysesp/reference/mar.md)
or
[`mnar()`](https://giovannitinervia9.github.io/baysesp/reference/mnar.md))
that states how the probability of being verified depends on the test
results and on the disease status. A model object holds the priors and
the verification mechanism; the data are supplied to
[`baysesp()`](https://giovannitinervia9.github.io/baysesp/reference/baysesp.md),
which fits the model.

## See also

[`saturated()`](https://giovannitinervia9.github.io/baysesp/reference/saturated.md),
[`loglinear()`](https://giovannitinervia9.github.io/baysesp/reference/loglinear.md),
[`random_effects()`](https://giovannitinervia9.github.io/baysesp/reference/random_effects.md)
