# Verification mechanisms

`verification_mechanism` is the abstract parent class of the assumptions
on how subjects are selected for verification by the reference standard.
It cannot be instantiated: the concrete assumptions are the missing at
random assumption, in which verification depends only on the test
results, and the missing not at random assumption, in which verification
also depends on the disease status (Rubin, 1976; Little and Rubin,
2002). Each concrete assumption has its own page.

## Usage

``` r
verification_mechanism()
```

## References

Little, R. J. A. and Rubin, D. B. (2002). *Statistical Analysis with
Missing Data* (2nd ed.). Wiley.

Rubin, D. B. (1976). Inference and missing data. *Biometrika*, 63(3),
581–592.
[doi:10.1093/biomet/63.3.581](https://doi.org/10.1093/biomet/63.3.581)

## See also

[`mar()`](https://giovannitinervia9.github.io/baysesp/reference/mar.md),
[`mnar()`](https://giovannitinervia9.github.io/baysesp/reference/mnar.md)

## Examples

``` r
S7::S7_inherits(mar(), verification_mechanism)
#> [1] TRUE
```
