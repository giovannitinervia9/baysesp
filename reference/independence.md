# Independent tests in simulated data

`independence()` is the assumption, used when simulating data, that the
test results are independent given the disease status: among the
diseased, and separately among the healthy, the result of one test
carries no information on the results of the others. The function has no
arguments. The object is passed to
[`simulate_study()`](https://giovannitinervia9.github.io/baysesp/reference/simulate_study.md)
as its `dependence` argument, where it is the default.

## Usage

``` r
independence()
```

## Value

An object of class `independence`.

## Details

Let \\t = (t_1, \dots, t_K)\\ be a response profile, with \\t_k = 1\\
for a positive result, and let \\Se_k\\ and \\Sp_k\\ be the sensitivity
and specificity of test \\k\\. The probabilities of the profiles among
the diseased and among the healthy are the products \$\$P(T = t \mid D
= 1) = \prod\_{k=1}^{K} Se_k^{t_k} (1 - Se_k)^{1 - t_k}, \qquad P(T = t
\mid D = 0) = \prod\_{k=1}^{K} (1 - Sp_k)^{t_k} Sp_k^{1 - t_k}.\$\$

## See also

[`simulate_study()`](https://giovannitinervia9.github.io/baysesp/reference/simulate_study.md),
[`latent_dependence()`](https://giovannitinervia9.github.io/baysesp/reference/latent_dependence.md),
[`probit_dependence()`](https://giovannitinervia9.github.io/baysesp/reference/probit_dependence.md)

## Examples

``` r
independence()
#> <independence> tests independent given the disease status
d <- simulate_study(n = 100, prevalence = 0.2, sensitivity = c(0.9, 0.8),
                    specificity = c(0.9, 0.95), verification_prob = 1,
                    dependence = independence())
```
