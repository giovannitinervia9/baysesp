# Random effects model

Models the joint distribution of the test results with one continuous
latent variable \\Z \sim N(0, 1)\\, independent of the disease and that
can be read as the severity of the condition. Within each disease class
the tests are independent given \\Z\\, and they are dependent overall
because they share \\Z\\. The model needs at least two tests.

## Usage

``` r
random_effects(
  prevalence_prior = prior_beta(0.5, 2),
  sensitivity_prior = prior_beta(0.5, 2),
  specificity_prior = prior_beta(0.5, 2),
  loadings_prior = prior_normal(0, 1),
  positive_loadings = TRUE,
  verification = mar()
)
```

## Arguments

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
  distribution on (0, 1). They are the priors of the marginal
  sensitivity and specificity, see Details.

- loadings_prior:

  Prior on the loadings \\b\_{dk}\\, a real-valued prior applied to all
  \\2K\\ loadings with the same prior, for example
  [`prior_normal()`](https://giovannitinervia9.github.io/baysesp/reference/prior_normal.md).
  If `positive_loadings = TRUE` the prior is restricted to positive
  values for all loadings; if `FALSE` it is restricted to positive
  values only for the first loading of each class, the others are real.
  Default `prior_normal(0, 1)`.

- positive_loadings:

  A single logical value, `TRUE` or `FALSE`, not `NA`. If `TRUE`
  (default), all the loadings are positive, so the tests can only be
  positively dependent within each class. If `FALSE`, only the first
  loading of each class is positive and the others can take any sign,
  which can make the sampler unreliable. See Details for why a
  constraint is needed and how the two options behave.

- verification:

  Verification mechanism: an object of class `mar` (default,
  [`mar()`](https://giovannitinervia9.github.io/baysesp/reference/mar.md))
  or
  [`mnar()`](https://giovannitinervia9.github.io/baysesp/reference/mnar.md).
  The priors of the verification mechanism (for example
  `verification = mnar(odds_ratio_prior = prior_normal(0, 1), by_profile = TRUE)`)
  are set in that object, not in this function.

## Value

An object of class `random_effects`.

## Details

Let \\t = (t_1, \dots, t_K)\\ be a response profile, with \\t_k = 1\\
for a positive result, and \\d \in \\0, 1\\\\ the disease status. Given
\\D = d\\ and \\Z = z\\, the test results are independent and \$\$P(T_k
= 1 \mid D = d, Z = z) = \Phi(a\_{dk} + b\_{dk} z),\$\$ where \\\Phi\\
is the standard normal distribution function, \\a\_{dk}\\ is an
intercept and \\b\_{dk}\\ a loading, for each test \\k\\ and class
\\d\\. The model has \\2K\\ intercepts and loadings per class, plus the
prevalence \\\pi = P(D = 1)\\. This is the Gaussian random effects
(probit latent variable) model for conditional dependence between tests
(Qu et al., 1996; Dendukuri and Joseph, 2001).

Integrating out \\Z\\ gives the probability of each profile in class
\\d\\, \\\theta\_{d,t} = P(T = t \mid D = d)\\: \$\$\theta\_{d,t} = \int
\prod\_{k=1}^{K} \Phi(a\_{dk} + b\_{dk} z)^{t_k}\\ \[1 - \Phi(a\_{dk} +
b\_{dk} z)\]^{1 - t_k}\\\phi(z)\\dz,\$\$ where \\\phi\\ is the standard
normal density. The integral is computed with 32-node Gauss-Hermite
quadrature.

Integrating out \\Z\\ also gives the marginal probits \$\$P(T_k = 1 \mid
D = d) = \Phi(m\_{dk}), \qquad m\_{dk} = \frac{a\_{dk}}{\sqrt{1 +
b\_{dk}^2}},\$\$ so that the sensitivity of test \\k\\ is \\Se_k =
\Phi(m\_{1k})\\ and its specificity is \\Sp_k = 1 - \Phi(m\_{0k})\\. The
model is sampled in the marginal probits \\m\_{dk}\\, which the data
inform directly, and \\a\_{dk} = m\_{dk}\sqrt{1 + b\_{dk}^2}\\. Since
\\Se_k = \Phi(m\_{1k})\\ and \\Sp_k = 1 - \Phi(m\_{0k})\\ whatever the
loadings, `sensitivity_prior` and `specificity_prior` are translated
into priors on \\m\_{1k}\\ and \\m\_{0k}\\ and are exactly the priors of
the marginal sensitivity and specificity. The uniform default
corresponds to a standard normal prior on \\m\_{dk}\\.

The loadings describe the dependence between the tests within each
class. Given \\D = d\\, the covariance between tests \\j\\ and \\k\\ has
the sign of \\b\_{dj} b\_{dk}\\: a positive product means the tests tend
to agree beyond what the disease explains, a negative product means they
tend to disagree. The size of the loadings sets the strength of the
dependence, and a loading close to 0 makes its test almost independent
of the others in that class. Loadings are specific to each class, so the
dependence can differ between diseased and healthy subjects.

Let \\n\_{t,1}\\ and \\n\_{t,0}\\ be the numbers of verified diseased
and verified healthy subjects with profile \\t\\, \\u_t\\ the number of
unverified subjects with profile \\t\\, \\\pi_1 = \pi\\ and \\\pi_0 =
1 - \pi\\. The log-likelihood of the aggregated data is \$\$\ell =
\sum_t \Big\[\sum\_{d} n\_{t,d} \log(\pi_d\\\theta\_{d,t}\\
\lambda\_{t,d}) + u_t \log \sum\_{d} \pi_d\\\theta\_{d,t}\\(1 -
\lambda\_{t,d})\Big\],\$\$ where \\\lambda\_{t,d} = P(V = 1 \mid T = t,
D = d)\\ is the probability of verification. Dependence models with
partial verification of this kind are discussed by Albert and Dodd
(2008).

With `verification = mar()` (missing at random) the verification
probability does not depend on the disease status (Rubin, 1976), so the
verification terms carry no information on the parameters and are
dropped: \$\$\ell = \sum_t \Big\[\sum\_{d} n\_{t,d}
\log(\pi_d\\\theta\_{d,t}) + u_t \log \sum\_{d}
\pi_d\\\theta\_{d,t}\Big\] + \textrm{const.}\$\$

With `verification = mnar()` (missing not at random) the verification
probability depends on the disease status (Kosinski and Barnhart, 2003):
\\\mathrm{logit}\\\lambda\_{t,d} = \alpha_t + \beta d\\ for every
profile with at least one verified subject, while \\\lambda\_{t,d} = 0\\
for profiles that are never verified. The parameter \\\beta\\ is the log
odds ratio of verification between diseased and healthy subjects with
the same profile; \\\beta \> 0\\ means that the diseased are verified
more often, and \\\exp(-\beta)\\ is the odds ratio of disease between
unverified and verified subjects with the same profile. With
`by_profile = TRUE` in `verification = mnar(by_profile = TRUE)`, each
profile has its own \\\beta_t\\.

The data do not inform every parameter directly. The probabilities
\\\theta\_{d,t}\\ of the profiles that are never verified are determined
only by the assumed latent structure, and under MNAR the parameter
\\\beta\\ is also informed only by that structure. Both depend on the
correctness of the model, so the results on these quantities are fragile
if the structure is wrong (Albert and Dodd, 2004).

**Constraint on the signs of the loadings.** Since \\Z\\ and \\-Z\\ have
the same distribution, changing the sign of all the loadings of a class
gives exactly the same profile probabilities: the parameters \\(a\_{d1},
\dots, a\_{dK}, b\_{d1}, \dots, b\_{dK})\\ and \\(a\_{d1}, \dots,
a\_{dK}, -b\_{d1}, \dots, -b\_{dK})\\ fit the data equally well. Without
a constraint the posterior would have two mirror-image modes, and
different chains could settle in different ones. Two constraints are
available:

- `positive_loadings = TRUE` (default): all the loadings are positive.
  The mirror image is excluded entirely, and every pair of tests is
  positively dependent within each class: the tests tend to err
  together, as is usual for tests that measure the same condition.
  Negative dependence cannot be represented, and independence (all
  loadings 0) is on the boundary of the parameter space: when the tests
  are in fact conditionally independent, the posterior keeps some
  positive dependence. The accuracy of each test is then estimated well,
  but quantities driven by the dependence are not: in simulations with
  independent tests, the 95% interval for the sensitivity of the rule
  "any test positive" covered the true value in only about 75% of the
  datasets, because the true value lies at the edge of what the model
  can represent. Compare the model with
  `loglinear(interaction_order = 1)` in
  [`compare_models()`](https://giovannitinervia9.github.io/baysesp/reference/compare_models.md)
  before relying on such quantities.

- `positive_loadings = FALSE`: only the first loading of each class is
  positive, the others can take any sign, so any sign of dependence is
  allowed. This removes the mirror image only when the first loading is
  clearly away from 0. When it is close to 0, for example because the
  first test is nearly independent of the others, the mirror-image
  configuration (first loading near 0, the others with the opposite
  sign) can be reached continuously, and the posterior keeps a second
  mode. Chains that settle in different modes give different estimates
  (R-hat well above 1), and the sampler can produce divergent
  transitions near the boundary. On simulated data with dependent tests
  this happens in most fits. Use this option only when negative
  dependence is plausible, with several chains, and check that the
  chains agree.

## References

Albert, P. S. and Dodd, L. E. (2004). A cautionary note on the
robustness of latent class models for estimating diagnostic error
without a gold standard. *Biometrics*, 60(2), 427–435.
[doi:10.1111/j.0006-341X.2004.00187.x](https://doi.org/10.1111/j.0006-341X.2004.00187.x)

Albert, P. S. and Dodd, L. E. (2008). On estimating diagnostic accuracy
from studies with multiple raters and partial gold standard evaluation.
*Journal of the American Statistical Association*, 103(481), 61–73.
[doi:10.1198/016214507000000329](https://doi.org/10.1198/016214507000000329)

Dendukuri, N. and Joseph, L. (2001). Bayesian approaches to modeling the
conditional dependence between multiple diagnostic tests. *Biometrics*,
57(1), 158–167.
[doi:10.1111/j.0006-341X.2001.00158.x](https://doi.org/10.1111/j.0006-341X.2001.00158.x)

Kosinski, A. S. and Barnhart, H. X. (2003). Accounting for nonignorable
verification bias in assessment of diagnostic tests. *Biometrics*,
59(1), 163–171.
[doi:10.1111/1541-0420.00019](https://doi.org/10.1111/1541-0420.00019)

Qu, Y., Tan, M. and Kutner, M. H. (1996). Random effects models in
latent class analysis for evaluating accuracy of diagnostic tests.
*Biometrics*, 52(3), 797–810.
[doi:10.2307/2533043](https://doi.org/10.2307/2533043)

Rubin, D. B. (1976). Inference and missing data. *Biometrika*, 63(3),
581–592.
[doi:10.1093/biomet/63.3.581](https://doi.org/10.1093/biomet/63.3.581)

## Examples

``` r
random_effects()
#> <random_effects> random effects model, Gauss-Hermite quadrature
#>   prevalence_prior:  Beta(mean = 0.5, precision = 2)
#>   sensitivity_prior: Beta(mean = 0.5, precision = 2)
#>   specificity_prior: Beta(mean = 0.5, precision = 2)
#>   loadings_prior:    Normal(mean = 0, sd = 1), all loadings positive
#>   verification:      mar
random_effects(positive_loadings = FALSE)
#> <random_effects> random effects model, Gauss-Hermite quadrature
#>   prevalence_prior:  Beta(mean = 0.5, precision = 2)
#>   sensitivity_prior: Beta(mean = 0.5, precision = 2)
#>   specificity_prior: Beta(mean = 0.5, precision = 2)
#>   loadings_prior:    Normal(mean = 0, sd = 1), first loading positive
#>   verification:      mar
```
