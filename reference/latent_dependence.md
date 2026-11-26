# Tests dependent through a shared latent variable in simulated data

`latent_dependence()` makes the simulated test results dependent given
the disease status through one continuous latent variable shared by all
the tests of a subject, such as the severity of the condition among the
diseased or a condition that causes false positives among the healthy.
The sensitivity and specificity of each test do not change: the latent
variable changes only how often the tests are right or wrong together.
The object is passed to
[`simulate_study()`](https://giovannitinervia9.github.io/baysesp/reference/simulate_study.md)
as its `dependence` argument.

## Usage

``` r
latent_dependence(diseased_loadings = 0, healthy_loadings = 0)
```

## Arguments

- diseased_loadings:

  Loadings \\b\_{1k}\\ of the tests among the diseased: a numeric vector
  of finite values, of any sign, either of length 1 (the same loading
  for all tests) or of length \\K\\ (one loading per test, in the order
  of the `sensitivity` argument of
  [`simulate_study()`](https://giovannitinervia9.github.io/baysesp/reference/simulate_study.md)).
  Any combination of finite values is allowed; the Details describe the
  dependences each combination produces. Default `0`: the tests are
  independent among the diseased. The length is checked against \\K\\
  when the data are simulated.

- healthy_loadings:

  Loadings \\b\_{0k}\\ of the tests among the healthy, with the same
  form as `diseased_loadings`. Default `0`: the tests are independent
  among the healthy.

## Value

An object of class `latent_dependence`.

## Details

Each subject has a latent variable \\Z \sim N(0, 1)\\, independent of
the disease status \\D\\. Given \\D = d\\ and \\Z = z\\ the \\K\\ tests
are independent, with \$\$P(T_k = 1 \mid D = d, Z = z) = \Phi(a\_{dk} +
b\_{dk} z),\$\$ where \\\Phi\\ is the standard normal distribution
function and \\b\_{dk}\\ is the loading of test \\k\\ in class \\d\\
(`diseased_loadings` for \\d = 1\\, `healthy_loadings` for \\d = 0\\).
The intercepts are \$\$a\_{1k} = \Phi^{-1}(Se_k)\sqrt{1 + b\_{1k}^2},
\qquad a\_{0k} = \Phi^{-1}(1 - Sp_k)\sqrt{1 + b\_{0k}^2},\$\$ so that,
integrating out \\Z\\, \\P(T_k = 1 \mid D = 1) = Se_k\\ and \\P(T_k = 0
\mid D = 0) = Sp_k\\ exactly, for any value of the loadings.

Equivalently, test \\k\\ is positive when a normal latent variable \\Y_k
= b\_{dk} Z + \varepsilon_k\\, with \\\varepsilon_k \sim N(0, 1)\\
independent, exceeds a threshold. Within class \\d\\ the correlation of
\\Y_j\\ and \\Y_k\\ is \$\$\rho\_{d,jk} = \frac{b\_{dj}\\
b\_{dk}}{\sqrt{(1 + b\_{dj}^2)(1 + b\_{dk}^2)}}.\$\$ Hence:

- a loading 0 makes its test independent of all the others in that
  class;

- tests \\j\\ and \\k\\ agree more often than under independence when
  \\b\_{dj} b\_{dk} \> 0\\, and less often when \\b\_{dj} b\_{dk} \<
  0\\;

- larger loadings in absolute value give stronger dependence. For
  example, two tests with sensitivity 0.8 and loadings 1 among the
  diseased have \\\rho = 0.5\\ and an odds ratio of about 4.7 between
  their results; with loadings 2, \\\rho = 0.8\\ and an odds ratio of
  about 19.

**Allowed values.** Every combination of finite real values is allowed:
the loadings of a class are not constrained by each other, by the
loadings of the other class or by the sensitivities and specificities,
and every combination gives a valid distribution of the results with the
requested sensitivities and specificities. The two classes are
independent choices, so everything below holds separately among the
diseased and among the healthy.

**Dependences that can be obtained.** Since all the tests share one
latent variable, the loadings cannot produce every pattern of pairwise
dependences. Write \\s_k = \|b\_{dk}\| / \sqrt{1 + b\_{dk}^2}\\, a
number in \\\[0, 1)\\ that measures how strongly test \\k\\ depends on
\\Z\\; then \\\|\rho\_{d,jk}\| = s_j s_k\\. Hence:

- **strength:** each test has a single strength \\s_k\\, and the
  dependence of a pair is the product of the strengths of its two tests.
  For three tests \\\|\rho\_{d,23}\| = \|\rho\_{d,12}\\ \rho\_{d,13}\| /
  s_1^2 \ge \|\rho\_{d,12}\\ \rho\_{d,13}\|\\, so tests 2 and 3 cannot
  be much less dependent on each other than both are on test 1: for
  example \\\rho\_{d,12} = \rho\_{d,13} = 0.8\\ forces \\\rho\_{d,23}
  \ge 0.64\\. All the correlations are strictly between \\-1\\ and 1;

- **sign:** the sign of \\\rho\_{d,jk}\\ is the product of the signs of
  \\b\_{dj}\\ and \\b\_{dk}\\. The tests split into those with positive
  and those with negative loading, with positive dependence within each
  group and negative dependence between the groups. Changing the sign of
  all the loadings of a class gives the same distribution, so only the
  relative signs matter. For any three tests the product
  \\\rho\_{d,12}\\ \rho\_{d,13}\\ \rho\_{d,23}\\ is never negative: for
  example, three tests cannot all be negatively dependent on each other;

- **independence:** \\\rho\_{d,jk} = 0\\ only when \\b\_{dj} = 0\\ or
  \\b\_{dk} = 0\\, and a loading 0 makes its test independent of all the
  others. Two tests cannot be independent of each other while both
  depend on a third.

To give the same correlation \\\rho \in \[0, 1)\\ to every pair, use the
same loading \\b = \sqrt{\rho / (1 - \rho)}\\ for all tests: \\b = 0.5\\
gives \\\rho = 0.2\\, \\b = 1\\ gives \\\rho = 0.5\\ and \\b = 2\\ gives
\\\rho = 0.8\\. Patterns outside these limits, such as a correlation of
any sign and size for each pair, are obtained with
`dependence = probit_dependence(...)`.

## See also

[`simulate_study()`](https://giovannitinervia9.github.io/baysesp/reference/simulate_study.md),
[`independence()`](https://giovannitinervia9.github.io/baysesp/reference/independence.md),
[`probit_dependence()`](https://giovannitinervia9.github.io/baysesp/reference/probit_dependence.md),
[`random_effects()`](https://giovannitinervia9.github.io/baysesp/reference/random_effects.md),
which fits the same latent structure

## Examples

``` r
# same loading for all tests among the diseased, independence among the healthy
latent_dependence(diseased_loadings = 1)
#> <latent_dependence> tests dependent through a shared latent variable
#>   diseased_loadings: 1
#>   healthy_loadings: 0

# one loading per test and class: A and B strongly dependent among the
# diseased, C almost independent; among the healthy only A and B, weakly
dep <- latent_dependence(diseased_loadings = c(1, 1, 0.3),
                         healthy_loadings = c(0.5, 0.5, 0))
d <- simulate_study(n = 2000, prevalence = 0.3,
                    sensitivity = c(A = 0.9, B = 0.8, C = 0.7),
                    specificity = c(A = 0.85, B = 0.9, C = 0.95),
                    verification_prob = 1, dependence = dep)
sick <- d$disease == 1
cor(d[sick, c("A", "B", "C")])
#>           A         B         C
#> A 1.0000000 0.2771135 0.1627437
#> B 0.2771135 1.0000000 0.1122374
#> C 0.1627437 0.1122374 1.0000000

# mixed signs: A and B agree, C disagrees with both
latent_dependence(diseased_loadings = c(1, 1, -1))
#> <latent_dependence> tests dependent through a shared latent variable
#>   diseased_loadings:  1,  1, -1
#>   healthy_loadings: 0
```
