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
  sensitivity_prior = prior_beta(0.5, 2),
  specificity_prior = prior_beta(0.5, 2),
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

- sensitivity_prior, specificity_prior:

  Priors on the sensitivity and on the specificity of the tests: a
  single prior on a probability,
  [`prior_beta()`](https://giovannitinervia9.github.io/baysesp/reference/prior_beta.md)
  or
  [`elicit_beta()`](https://giovannitinervia9.github.io/baysesp/reference/elicit_beta.md),
  used for all tests, or a list with one such prior per test, in the
  order of the tests. Default `prior_beta(0.5, 2)`, the uniform
  distribution on (0, 1). With `interaction_order` of 2 or more they act
  on conditional sensitivities and specificities: see the section
  "Priors on sensitivity and specificity".

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
normalisation runs over all \\2^K\\ profiles. Log-linear models of this
kind are used for multiple tests with partial verification (Baker,
1995).

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
0\\). Ignoring the dependence between tests given the disease status
biases the accuracy (Vacek, 1985; Torrance-Rynard and Walter, 1997;
Dendukuri and Joseph, 2001).

The log-likelihood of the data, with \\n\_{t,d}\\ verified subjects of
profile \\t\\ and disease status \\d\\, \\u_t\\ unverified subjects of
profile \\t\\, \\\pi_1 = \pi\\ the prevalence and \\\pi_0 = 1 - \pi\\,
is \$\$\ell = \sum_t \Big\[\sum_d n\_{t,d}
\log(\pi_d\\\theta\_{d,t}\\\lambda\_{t,d}) + u_t \log \sum_d
\pi_d\\\theta\_{d,t}\\(1 - \lambda\_{t,d})\Big\],\$\$ where
\\\lambda\_{t,d}\\ is the probability of verification. Dependence models
with partial verification of this kind are discussed by Albert and Dodd
(2008). With `verification = mar()` the verification probability does
not depend on the disease status given the profile (Rubin, 1976), so the
verification term is separated from the disease model and can be
ignored: \$\$\ell = \sum_t \Big\[\sum_d n\_{t,d}
\log(\pi_d\\\theta\_{d,t}) + u_t \log \sum_d
\pi_d\\\theta\_{d,t}\Big\]\$\$ plus a term that depends only on the
verification probabilities. With `verification = mnar()` the
verification probability may depend on the disease status (Baker, 1995;
Kosinski and Barnhart, 2003): \$\$\mathrm{logit}\\\lambda\_{t,d} =
\alpha_t + \beta d,\$\$ with one \\\alpha_t\\ for each verified profile.
The parameter \\\beta\\ is the log odds ratio of verification between
diseased and healthy subjects with the same profile, and
\\\exp(-\beta)\\ is the odds ratio of disease between unverified and
verified subjects with the same profile. With
`verification = mnar(by_profile = TRUE)` there is one \\\beta_t\\ for
each profile. Profiles that are never verified have \\\lambda\_{t,d} =
0\\.

With order lower than \\K\\ the restricted structure also determines the
disease proportion in the never-verified profiles, and makes \\\beta\\
formally estimable. This information comes from the assumed structure,
not from the data, and it is fragile if the structure is wrong (Albert
and Dodd, 2004). With order \\K\\ the model imposes no restriction and
\\\beta\\ is not identified (Gustafson, 2015).

## Priors on sensitivity and specificity

The model is fitted on the main effects \\\eta\_{d,k}\\, and
`sensitivity_prior` and `specificity_prior` are translated into priors
on them. With `interaction_order = 1` the main effects are the logits of
the accuracy, \$\$\eta\_{1,k} = \mathrm{logit}\\\mathrm{Se}\_k, \qquad
\eta\_{0,k} = \mathrm{logit}(1 - \mathrm{Sp}\_k),\$\$ so a prior
\\\mathrm{Beta}(a, b)\\ given to `sensitivity_prior` is the prior of
\\\mathrm{Se}\_k\\, and one given to `specificity_prior` is the prior of
\\\mathrm{Sp}\_k\\ (that is, \\\mathrm{Beta}(b, a)\\ on the false
positive rate \\1 - \mathrm{Sp}\_k\\).

With `interaction_order` of 2 or more the main effects are logits of
conditional probabilities: the probability that test \\k\\ is positive
among the subjects of class \\d\\ who are negative on all the other
tests, \$\$\mathrm{expit}\\\eta\_{d,k} = P(T_k = 1 \mid D = d,\\ T_j = 0
\\ (j \neq k)).\$\$ So `sensitivity_prior` acts on the conditional
sensitivity \$\$P(T_k = 1 \mid D = 1,\\ T_j = 0 \\ (j \neq k)) =
\mathrm{expit}\\\eta\_{1,k},\$\$ and `specificity_prior` on the
conditional specificity \$\$P(T_k = 0 \mid D = 0,\\ T_j = 0 \\ (j \neq
k)) = 1 - \mathrm{expit}\\\eta\_{0,k},\$\$ not on the marginal
\\\mathrm{Se}\_k\\ and \\\mathrm{Sp}\_k\\. The priors implied on the
marginal sensitivity and specificity also depend on
`interactions_prior`. With the default priors and three tests they are
close to the uniform, but the interactions move the marginal sensitivity
slightly up and the marginal specificity slightly down (medians about
0.54 and 0.46). An informative `sensitivity_prior` or
`specificity_prior` should be set with this in mind.

## References

Albert, P. S. and Dodd, L. E. (2004). A cautionary note on the
robustness of latent class models for estimating diagnostic error
without a gold standard. *Biometrics*, 60(2), 427–435.
[doi:10.1111/j.0006-341X.2004.00187.x](https://doi.org/10.1111/j.0006-341X.2004.00187.x)

Albert, P. S. and Dodd, L. E. (2008). On estimating diagnostic accuracy
from studies with multiple raters and partial gold standard evaluation.
*Journal of the American Statistical Association*, 103(481), 61–73.
[doi:10.1198/016214507000000329](https://doi.org/10.1198/016214507000000329)

Baker, S. G. (1995). Evaluating multiple diagnostic tests with partial
verification. *Biometrics*, 51(1), 330–337.
[doi:10.2307/2533339](https://doi.org/10.2307/2533339)

Dendukuri, N. and Joseph, L. (2001). Bayesian approaches to modeling the
conditional dependence between multiple diagnostic tests. *Biometrics*,
57(1), 158–167.
[doi:10.1111/j.0006-341X.2001.00158.x](https://doi.org/10.1111/j.0006-341X.2001.00158.x)

Gustafson, P. (2015). *Bayesian Inference for Partially Identified
Models: Exploring the Limits of Limited Data*. Chapman & Hall/CRC, Boca
Raton, FL.

Kosinski, A. S. and Barnhart, H. X. (2003). Accounting for nonignorable
verification bias in assessment of diagnostic tests. *Biometrics*,
59(1), 163–171.
[doi:10.1111/1541-0420.00019](https://doi.org/10.1111/1541-0420.00019)

Rubin, D. B. (1976). Inference and missing data. *Biometrika*, 63(3),
581–592.
[doi:10.1093/biomet/63.3.581](https://doi.org/10.1093/biomet/63.3.581)

Torrance-Rynard, V. L. and Walter, S. D. (1997). Effects of dependent
errors in the assessment of diagnostic test performance. *Statistics in
Medicine*, 16(19), 2157–2175.
[doi:10.1002/(SICI)1097-0258(19971015)16:19\<2157::AID-SIM653\>3.0.CO;2-X](https://doi.org/10.1002/%28SICI%291097-0258%2819971015%2916%3A19%3C2157%3A%3AAID-SIM653%3E3.0.CO%3B2-X)

Vacek, P. M. (1985). The effect of conditional dependence on the
evaluation of diagnostic tests. *Biometrics*, 41(4), 959–968.
[doi:10.2307/2530967](https://doi.org/10.2307/2530967)

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
#>   sensitivity_prior:  Beta(mean = 0.5, precision = 2)
#>   specificity_prior:  Beta(mean = 0.5, precision = 2)
#>   interactions_prior: Normal(mean = 0, sd = 1)
#>   verification:       mar
loglinear(interaction_order = 2)  # pairwise dependence
#> <loglinear> log-linear model of order 2
#>   prevalence_prior:   Beta(mean = 0.5, precision = 2)
#>   sensitivity_prior:  Beta(mean = 0.5, precision = 2)
#>   specificity_prior:  Beta(mean = 0.5, precision = 2)
#>   interactions_prior: Normal(mean = 0, sd = 1)
#>   verification:       mar
```
