# Beta prior with mean and precision

A Beta distribution for a probability, written in terms of its mean and
its precision (the sum of the two shape parameters). The usual shape
parameters are `shape1 = mean * precision` and
`shape2 = (1 - mean) * precision`; they are available as the read-only
properties `@shape1` and `@shape2`.

## Usage

``` r
prior_beta(mean = 0.5, precision = 2)
```

## Arguments

- mean:

  Mean of the distribution: a single number in (0, 1). Default 0.5.

- precision:

  Precision: a single positive number. The larger it is, the more
  concentrated the distribution around its mean. The default 2, with
  `mean = 0.5`, gives the uniform distribution on (0, 1).

## Value

An object of class `prior_beta`.

## See also

[`elicit_beta()`](https://giovannitinervia9.github.io/baysesp/reference/elicit_beta.md)
to obtain the precision from a plausible upper value.

## Examples

``` r
prior_beta(mean = 0.5, precision = 2)   # uniform
#> <prior_beta> Beta(mean = 0.5, precision = 2)
#>   acts on a probability
p <- prior_beta(mean = 0.1, precision = 50)
p@shape1
#> [1] 5
p@shape2
#> [1] 45
```
