# Beta prior elicited from a mean and a plausible upper value

Builds a
[`prior_beta()`](https://giovannitinervia9.github.io/baysesp/reference/prior_beta.md)
from two simple statements about a probability: its most plausible value
(`mean`) and a value that it is unlikely to exceed (`upper`), in the
sense that `P(X <= upper) = prob`. For example, "the proportion of
diseased among never-verified subjects is about 0.5% and almost
certainly below 1.5%" is `elicit_beta(mean = 0.005, upper = 0.015)`.

## Usage

``` r
elicit_beta(mean, upper, prob = 0.95)
```

## Arguments

- mean:

  Most plausible value: a single number in (0, 1). No default.

- upper:

  Plausible upper value: a single number greater than `mean` and lower
  than 1. No default.

- prob:

  Probability that the value is below `upper`: a single number in (0.5,
  1). Default 0.95.

## Value

An object of class `prior_beta`.

## Details

The probability `P(X <= upper)` of a Beta distribution with a given mean
is not monotone in the precision: there are usually two solutions, a
U-shaped Beta with its mass at the extremes and a regular one. The
function returns the regular one, the root to the right of the minimum.

## Examples

``` r
p <- elicit_beta(mean = 0.005, upper = 0.015)
p
#> <prior_beta> Beta(mean = 0.005, precision = 197.8)
#>   acts on a probability
stats::pbeta(0.015, p@shape1, p@shape2)
#> [1] 0.95
```
