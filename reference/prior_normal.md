# Normal prior

A normal distribution for an unconstrained parameter. When it is given
to a probability, it acts on its logit.

## Usage

``` r
prior_normal(mean = 0, sd = 1)
```

## Arguments

- mean:

  Mean: a single finite number. Default 0.

- sd:

  Standard deviation: a single positive number. Default 1.

## Value

An object of class `prior_normal`.

## Examples

``` r
prior_normal(0, 1.5)
#> <prior_normal> Normal(mean = 0, sd = 1.5)
#>   acts on an unconstrained parameter (on the logit of a probability)
```
