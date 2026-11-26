# Tests with correlated latent variables in simulated data

`probit_dependence()` makes the simulated test results dependent given
the disease status by giving each test a normal latent variable and
letting these variables be correlated, with one correlation for each
pair of tests and disease class. The correlations can have any sign,
separately for each pair. The sensitivity and specificity of each test
do not change: the correlations change only how often the tests are
right or wrong together. The object is passed to
[`simulate_study()`](https://giovannitinervia9.github.io/baysesp/reference/simulate_study.md)
as its `dependence` argument.

## Usage

``` r
probit_dependence(diseased_correlation = 0, healthy_correlation = 0)
```

## Arguments

- diseased_correlation:

  Correlations \\R_1\\ between the latent variables of the tests among
  the diseased: either a single number, the same correlation for every
  pair of tests, or a symmetric \\K \times K\\ matrix with 1 on the
  diagonal, whose rows and columns follow the order of the `sensitivity`
  argument of
  [`simulate_study()`](https://giovannitinervia9.github.io/baysesp/reference/simulate_study.md).
  Values in \\\[-1, 1\]\\. Default `0`: the tests are independent among
  the diseased. The size of the matrix and its positive definiteness
  (which excludes, for example, a single number equal to 1 with two or
  more tests) are checked when the data are simulated.

- healthy_correlation:

  Correlations \\R_0\\ between the latent variables of the tests among
  the healthy, with the same form as `diseased_correlation`. Default
  `0`: the tests are independent among the healthy.

## Value

An object of class `probit_dependence`.

## Details

Given the disease status \\D = d\\, the \\K\\ tests have latent
variables \\Y = (Y_1, \dots, Y_K)\\ jointly normal with mean 0, variance
1 and correlation matrix \\R_d\\ (`diseased_correlation` for \\d = 1\\,
`healthy_correlation` for \\d = 0\\). Test \\k\\ is positive when
\$\$Y_k \< \Phi^{-1}(Se_k) \textrm{ if } d = 1, \qquad Y_k \<
\Phi^{-1}(1 - Sp_k) \textrm{ if } d = 0,\$\$ where \\\Phi\\ is the
standard normal distribution function, so that \\P(T_k = 1 \mid D = 1) =
Se_k\\ and \\P(T_k = 0 \mid D = 0) = Sp_k\\ exactly, whatever the
correlations.

The element \\R\_{d,jk}\\ is the correlation between the latent
variables of tests \\j\\ and \\k\\ in class \\d\\ (known as the
tetrachoric correlation of the two results). Positive values make the
two tests agree more often than under independence, negative values less
often, and 0 makes their results uncorrelated. For example, two tests
with sensitivity 0.8 and correlation 0.5 among the diseased have an odds
ratio of about 4.7 between their results; with correlation 0.8, about
19. With all the correlations of a class equal to 0 the tests are
independent in that class.

A correlation matrix must be positive definite: not every set of
pairwise correlations is possible. For example, the same correlation
\\r\\ for all pairs of \\K\\ tests needs \\-1/(K - 1) \< r \< 1\\, so
with three tests \\r \> -1/2\\. Since \\K\\ is known only when the data
are simulated, this condition is checked by
[`simulate_study()`](https://giovannitinervia9.github.io/baysesp/reference/simulate_study.md),
which gives an error for a matrix that is not positive definite.

## See also

[`simulate_study()`](https://giovannitinervia9.github.io/baysesp/reference/simulate_study.md),
[`independence()`](https://giovannitinervia9.github.io/baysesp/reference/independence.md),
[`latent_dependence()`](https://giovannitinervia9.github.io/baysesp/reference/latent_dependence.md)

## Examples

``` r
# correlation 0.4 between every pair among the diseased, 0.2 among the healthy
probit_dependence(diseased_correlation = 0.4, healthy_correlation = 0.2)
#> <probit_dependence> tests with correlated latent variables
#>   diseased_correlation: 0.4
#>   healthy_correlation: 0.2

# different correlation for each pair: A and B positively correlated,
# C negatively correlated with both
R1 <- matrix(c( 1,    0.5, -0.3,
                0.5,  1,   -0.3,
               -0.3, -0.3,  1), 3, 3)
d <- simulate_study(n = 2000, prevalence = 0.3,
                    sensitivity = c(A = 0.9, B = 0.8, C = 0.7),
                    specificity = c(A = 0.85, B = 0.9, C = 0.95),
                    verification_prob = 1,
                    dependence = probit_dependence(diseased_correlation = R1))
sick <- d$disease == 1
cor(d[sick, c("A", "B", "C")])
#>             A           B           C
#> A  1.00000000  0.16611978 -0.09133949
#> B  0.16611978  1.00000000 -0.07316651
#> C -0.09133949 -0.07316651  1.00000000

# every pair negatively correlated: possible with three tests if r > -1/2
probit_dependence(diseased_correlation = -0.3)
#> <probit_dependence> tests with correlated latent variables
#>   diseased_correlation: -0.3
#>   healthy_correlation: 0
```
