# Verification mechanisms

`verification_mechanism` is the abstract parent class of the assumptions
on how subjects are selected for verification by the reference standard.
It cannot be instantiated: the concrete assumptions are the missing at
random assumption, in which verification depends only on the test
results, and the missing not at random assumption, in which verification
also depends on the disease status. Each concrete assumption has its own
page.

## Usage

``` r
verification_mechanism()
```

## See also

[`mar()`](https://giovannitinervia9.github.io/baysesp/reference/mar.md),
[`mnar()`](https://giovannitinervia9.github.io/baysesp/reference/mnar.md)

## Examples

``` r
S7::S7_inherits(mar(), verification_mechanism)
#> [1] TRUE
```
