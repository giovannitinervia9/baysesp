
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

y <- baysesp_simul(n = 1000, p = 0.01, se = 0.9, sp = 0.5)
y
#>    nsim1
#> s1     7
#> s2     1
#> r1   149
#> r2   258
#> u1   355
#> u2   230
```

Next, we convert the simulated data into the appropriate format for the
`baysesp()` function using `y_to_stan_data()`:

``` r
stan_data <- y_to_stan_data(y)
stan_data
#> $s1
#> [1] 7
#> 
#> $s2
#> [1] 1
#> 
#> $r1
#> [1] 149
#> 
#> $r2
#> [1] 258
#> 
#> $u1
#> [1] 355
#> 
#> $u2
#> [1] 230
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
#>           mean se_mean    sd      2.5%     97.5%     n_eff Rhat
#> p        0.013   0.000 0.005     0.005     0.025 219100.67    1
#> se       0.745   0.000 0.150     0.400     0.964 246885.74    1
#> sp       0.492   0.000 0.016     0.461     0.524 313749.55    1
#> l11      0.704   0.000 0.095     0.503     0.869 292434.20    1
#> l21      0.312   0.000 0.020     0.273     0.352 320614.10    1
#> l12      0.475   0.000 0.156     0.185     0.780 267210.34    1
#> l22      0.527   0.000 0.023     0.483     0.571 312806.13    1
#> lp__ -1431.748   0.006 1.920 -1436.350 -1429.016  88972.24    1
```
