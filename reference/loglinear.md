# Log-linear model of the test profiles

Defines a model in which the distribution of the response profiles is
log-linear in each disease class, so the dependence between the tests is
restricted to interactions of at most `interaction_order` tests. The
object is passed to
[`baysesp()`](https://giovannitinervia9.github.io/baysesp/reference/baysesp.md)
as its `model` argument.

## Usage

``` r
loglinear(
  interaction_order = 1L,
  prevalence_prior = prior_beta(0.5, 2),
  main_effects_prior = prior_normal(0, 1.5),
  interactions_prior = prior_normal(0, 1),
  verification = mar()
)
```

## Arguments

- interaction_order:

  A single whole number from 1 to the number of tests \\K\\, the largest
  number of tests in an interaction. Default `1L` (conditional
  independence). The lower bound is checked when the model is created;
  the upper bound is checked when the model is fitted.

- prevalence_prior:

  Prior on the prevalence \\\pi\\: a single prior, either on a
  probability, such as
  [`prior_beta()`](https://giovannitinervia9.github.io/baysesp/reference/prior_beta.md),
  which acts on \\\pi\\, or on the real line, such as
  [`prior_normal()`](https://giovannitinervia9.github.io/baysesp/reference/prior_normal.md),
  which acts on \\\mathrm{logit}\\\pi\\. Default `prior_beta(0.5, 2)`,
  the uniform distribution on (0, 1).

- main_effects_prior:

  A single prior on the real line, applied to each class-specific main
  effect \\\eta\_{d,k}\\ (one per test and disease class). Default
  `prior_normal(0, 1.5)`. With order 1 these parameters are logits of
  sensitivities and of false positive rates.

- interactions_prior:

  A single prior on the real line, applied to each class-specific
  interaction of order 2 or higher. Default `prior_normal(0, 1)`. It has
  no effect when `interaction_order = 1`.

- verification:

  The verification mechanism:
  [`mar()`](https://giovannitinervia9.github.io/baysesp/reference/mar.md)
  (default), with verification that depends only on the test profile, or
  [`mnar()`](https://giovannitinervia9.github.io/baysesp/reference/mnar.md),
  with verification that may also depend on the disease status.

## Value

An object of class `loglinear`.

## Details

Let \\K\\ be the number of tests, \\t = (t_1, \dots, t_K)\\ a response
profile with \\t_k = 1\\ for a positive result, \\d = 1\\ for diseased
and \\d = 0\\ for healthy, and let \\m\\ be `interaction_order`. For
each disease class the profile probabilities \\\theta_d\\ are
\$\$\theta_d = \mathrm{softmax}(X_m \eta_d),\$\$ where \\X_m\\ is the
design matrix whose columns are the functions \\x_S(t) = \prod\_{k \in
S} t_k\\, one for each set \\S\\ of at most \\m\\ tests, and \\\eta_d\\
is the vector of the corresponding class-specific parameters
\\\eta\_{d,S}\\. Equivalently, \$\$\theta\_{d,t} \propto
\exp\Big(\sum\_{S} \eta\_{d,S}\\ x_S(t)\Big),\$\$ where the
normalisation runs over all \\2^K\\ profiles.

The order controls the structure of the dependence:

- Order 1 (`interaction_order = 1`): the tests are conditionally
  independent given the disease status.

- Order 2: each pair of tests may be associated within each class, but
  there are no interactions among three or more tests.

- Order \\K\\: no restriction, every profile distribution is allowed in
  each class.

The parameters have a direct interpretation. With order 1, \\\eta\_{1,k}
= \mathrm{logit}(Se_k)\\ and \\\eta\_{0,k} = \mathrm{logit}(1 - Sp_k)\\,
so the main effects are the logits of the sensitivity and of the false
positive rate of test \\k\\. A pairwise term \\\eta\_{d,\\j,k\\}\\ is
the log odds ratio between tests \\j\\ and \\k\\ within class \\d\\,
given the results of the other tests (with order 2 it is the same
whatever those results; with higher orders it is the value when the
other tests are negative): zero means that the two tests are independent
within that class, positive values mean that they are positively
associated among the diseased (\\d = 1\\) or among the healthy (\\d =
0\\).

The log-likelihood of the data, with \\n\_{t,d}\\ verified subjects of
profile \\t\\ and disease status \\d\\, \\u_t\\ unverified subjects of
profile \\t\\, \\\pi_1 = \pi\\ the prevalence and \\\pi_0 = 1 - \pi\\,
is \$\$\ell = \sum_t \Big\[\sum_d n\_{t,d}
\log(\pi_d\\\theta\_{d,t}\\\lambda\_{t,d}) + u_t \log \sum_d
\pi_d\\\theta\_{d,t}\\(1 - \lambda\_{t,d})\Big\],\$\$ where
\\\lambda\_{t,d}\\ is the probability of verification. With
`verification = mar()` the verification probability does not depend on
the disease status given the profile, so the verification term is
separated from the disease model and can be ignored: \$\$\ell = \sum_t
\Big\[\sum_d n\_{t,d} \log(\pi_d\\\theta\_{d,t}) + u_t \log \sum_d
\pi_d\\\theta\_{d,t}\Big\]\$\$ plus a term that depends only on the
verification probabilities. With `verification = mnar()` the
verification probability may depend on the disease status:
\$\$\mathrm{logit}\\\lambda\_{t,d} = \alpha_t + \beta d,\$\$ with one
\\\alpha_t\\ for each verified profile. The parameter \\\beta\\ is the
log odds ratio of verification between diseased and healthy subjects
with the same profile, and \\\exp(-\beta)\\ is the odds ratio of disease
between unverified and verified subjects with the same profile. With
`verification = mnar(by_profile = TRUE)` there is one \\\beta_t\\ for
each profile. Profiles that are never verified have \\\lambda\_{t,d} =
0\\.

With order lower than \\K\\ the restricted structure also determines the
disease proportion in the never-verified profiles, and makes \\\beta\\
formally estimable. This information comes from the assumed structure,
not from the data, and it is fragile if the structure is wrong. With
order \\K\\ the model imposes no restriction and \\\beta\\ is not
identified.

## See also

[`saturated()`](https://giovannitinervia9.github.io/baysesp/reference/saturated.md)
for the model without restrictions on the profile distribution,
[`mar()`](https://giovannitinervia9.github.io/baysesp/reference/mar.md)
and
[`mnar()`](https://giovannitinervia9.github.io/baysesp/reference/mnar.md)
for the verification mechanisms.

## Examples

``` r
loglinear()                       # conditional independence
#> <loglinear> log-linear model of order 1
#>   prevalence_prior:   Beta(mean = 0.5, precision = 2)
#>   main_effects_prior: Normal(mean = 0, sd = 1.5)
#>   interactions_prior: Normal(mean = 0, sd = 1)
#>   verification:       mar
loglinear(interaction_order = 2)  # pairwise dependence
#> <loglinear> log-linear model of order 2
#>   prevalence_prior:   Beta(mean = 0.5, precision = 2)
#>   main_effects_prior: Normal(mean = 0, sd = 1.5)
#>   interactions_prior: Normal(mean = 0, sd = 1)
#>   verification:       mar
```
