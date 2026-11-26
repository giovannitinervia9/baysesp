# Student t prior

A Student t distribution for an unconstrained parameter, with heavier
tails than the normal. When it is given to a probability, it acts on its
logit.

## Usage

``` r
prior_student(df = 3, location = 0, scale = 1)
```

## Arguments

- df:

  Degrees of freedom: a single positive number. Default 3.

- location:

  Location, the centre of the distribution: a single finite number.
  Default 0.

- scale:

  Scale: a single positive number. Default 1.

## Value

An object of class `prior_student`.

## Examples

``` r
prior_student(df = 3, location = 0, scale = 1)
#> <prior_student> Student-t(df = 3, location = 0, scale = 1)
#>   acts on an unconstrained parameter (on the logit of a probability)
```
