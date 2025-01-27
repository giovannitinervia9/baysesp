
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `baysesp`

<!-- badges: start -->
<!-- badges: end -->

`baysesp` is an `R` package that enables users to estimate the
sensitivity and specificity of diagnostic tests while accounting for
verification bias using Bayesian methods.

## Installation

You can install the development version of `baysesp` from
[GitHub](https://github.com/giovannitinervia9/baysesp) with:

``` r
devtools::install_github("giovannitinervia9/baysesp")
```

## Example

We first generate a dataset using the `baysesp_simul()` function:

``` r
library(baysesp)
set.seed(123)
y <- baysesp_simul(n = 1000, p = 0.01, se = 0.9, sp = 0.5)
```

Next, we convert the simulated data into the appropriate format for the
`baysesp()` function using `y_to_stan_data()`:

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
s <- summary(fit, probs = c(a / 2, 1 - a / 2))
round(s$summary, 3)
#>           mean se_mean    sd      2.5%     97.5%    n_eff  Rhat
#> p        0.012   0.000 0.005     0.004     0.024 4529.181 1.000
#> se       0.718   0.002 0.165     0.351     0.960 4991.095 0.999
#> sp       0.510   0.000 0.016     0.479     0.541 5827.087 1.000
#> l11      0.702   0.001 0.095     0.502     0.869 6253.651 1.000
#> l21      0.312   0.000 0.021     0.272     0.355 6553.005 1.000
#> l12      0.474   0.002 0.160     0.180     0.787 4830.410 1.000
#> l22      0.487   0.000 0.022     0.443     0.530 6432.206 1.000
#> lp__ -1429.970   0.045 1.932 -1434.735 -1427.176 1879.181 1.001
```
