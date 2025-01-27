
<!-- README.md is generated from README.Rmd. Please edit that file -->

# baysesp

<!-- badges: start -->
<!-- badges: end -->

`baysesp` is an `R` package that allows users to estimate the
sensitivity and specificity of diagnostic tests in the presence of
verification bias via bayesian methods.

## Installation

You can install the development version of `baysesp` from
[GitHub](https://github.com/giovannitinervia9/baysesp) with:

``` r
devtools::install_github("giovannitinervia9/baysesp")
```

## Example

We first generate a dataset using the `baysesp_simul()`

``` r
library(baysesp)
set.seed(123)
y <- baysesp_simul(n = 1000, p = 0.01, se = 0.9, sp = 0.5)
```

Next, we convert the simulated data into the appropriate format for the
`baysesp()` function using `y_to_stan_data()`

``` r
stan_data <- y_to_stan_data(y)
```

We specify our priors as follows:

``` r
priors <- list(
  p = "beta(0.0199, 1.9701)",
  se = "beta(1, 1)",
  sp = "beta(1, 1)",
  l11 = "beta(15.4, 6.6)",
  l21 = "beta(15, 10)",
  l12 = "beta(4.67, 4.67)",
  l22 = "beta(0.4, 3.6)"
)
```

We now fit the model using `baysesp()`:

``` r
fit <- baysesp(stan_data, priors = priors)
```

Let’s check the results:

``` r
a <- 0.05
s <- summary(fit)
round(s$summary, 3)
#>           mean se_mean    sd      2.5%       25%       50%       75%     97.5%
#> p        0.012   0.000 0.005     0.004     0.008     0.011     0.014     0.024
#> se       0.718   0.002 0.165     0.351     0.614     0.742     0.846     0.960
#> sp       0.510   0.000 0.016     0.479     0.499     0.509     0.521     0.541
#> l11      0.702   0.001 0.095     0.502     0.640     0.709     0.772     0.869
#> l21      0.312   0.000 0.021     0.272     0.298     0.312     0.326     0.355
#> l12      0.474   0.002 0.160     0.180     0.357     0.470     0.589     0.787
#> l22      0.487   0.000 0.022     0.443     0.471     0.486     0.502     0.530
#> lp__ -1429.970   0.045 1.932 -1434.735 -1431.123 -1429.659 -1428.531 -1427.176
#>         n_eff  Rhat
#> p    4529.181 1.000
#> se   4991.095 0.999
#> sp   5827.087 1.000
#> l11  6253.651 1.000
#> l21  6553.005 1.000
#> l12  4830.410 1.000
#> l22  6432.206 1.000
#> lp__ 1879.181 1.001
```
