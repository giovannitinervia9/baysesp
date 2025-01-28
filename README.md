
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

y <- baysesp_simul(n = 10000, p = 0.05, se = 0.8, sp = 0.9)
y
#>    nsim1
#> s1   322
#> s2    39
#> r1   271
#> r2  4269
#> u1   717
#> u2  4382
```

Next, we convert the simulated data into the appropriate format for the
`baysesp()` function using `y_to_stan_data()`:

``` r
stan_data <- y_to_stan_data(y)
stan_data
#> $s1
#> [1] 322
#> 
#> $s2
#> [1] 39
#> 
#> $r1
#> [1] 271
#> 
#> $r2
#> [1] 4269
#> 
#> $u1
#> [1] 717
#> 
#> $u2
#> [1] 4382
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
fit <- baysesp(stan_data, priors = priors, chains = 8, iter = 50000, cores = 8)
```

Let’s check the results:

``` r
a <- 0.05
s <- summary(fit, probs = c(a / 2, 1 - a / 2))
round(s$summary, 3)
#>            mean se_mean    sd       2.5%      97.5%    n_eff Rhat
#> p         0.068   0.000 0.014      0.047      0.098 45786.28    1
#> se        0.852   0.000 0.062      0.695      0.934 57338.37    1
#> sp        0.922   0.000 0.013      0.902      0.950 45623.67    1
#> l11       0.587   0.001 0.119      0.381      0.816 49817.85    1
#> l21       0.390   0.000 0.079      0.294      0.598 43438.67    1
#> l12       0.463   0.001 0.157      0.175      0.773 76119.03    1
#> l22       0.497   0.000 0.006      0.485      0.509 88340.80    1
#> lp__ -11490.767   0.007 1.878 -11495.290 -11488.092 72012.69    1
```
