# Missing not at random verification

`mnar()` is the assumption that verification depends on the test results
and also on the disease status: among subjects with the same profile,
diseased and healthy subjects can be verified with different
probabilities. Its arguments are the prior on the log odds ratio of
verification, the prior on the verification intercepts and the choice of
whether the odds ratio varies across profiles.

## Usage

``` r
mnar(
  odds_ratio_prior = prior_normal(0, 1),
  intercepts_prior = prior_normal(0, 2.5),
  by_profile = FALSE
)
```

## Arguments

- odds_ratio_prior:

  Prior on \\\beta\\, the log odds ratio of verification between
  diseased and healthy subjects with the same profile. It must act on
  the real line, such as
  [`prior_normal()`](https://giovannitinervia9.github.io/baysesp/reference/prior_normal.md)
  or
  [`prior_student()`](https://giovannitinervia9.github.io/baysesp/reference/prior_student.md);
  a prior on a probability, such as
  [`prior_beta()`](https://giovannitinervia9.github.io/baysesp/reference/prior_beta.md),
  is not accepted. The default `prior_normal(0, 1)` centres \\\beta\\ at
  0, which is missing at random, and gives 95% of its prior mass to odds
  ratios between about 0.14 and 7.1.

- intercepts_prior:

  Prior on the verification intercepts \\\alpha_t\\, one for each
  profile with at least one verified subject. It must act on the real
  line, as for `odds_ratio_prior`. It is used only by the log-linear and
  the random-effects models; the saturated model does not use it. The
  default is `prior_normal(0, 2.5)`.

- by_profile:

  A single logical value, `TRUE` or `FALSE`, with default `FALSE`. If
  `FALSE`, one \\\beta\\ is shared by all profiles with at least one
  verified subject, so the odds ratio of verification is the same for
  every profile. If `TRUE`, each of these profiles has its own
  \\\beta_t\\, all with prior `odds_ratio_prior`, so the odds ratio may
  differ across profiles. Profiles with no verified subject have no
  \\\beta\\.

## Value

An object of class `mnar`.

## Details

Let \\V = 1\\ indicate that a subject is verified, \\T\\ be the response
profile, \\D\\ the disease status, and \\\lambda\_{t,d} = P(V = 1 \mid T
= t, D = d)\\. The model is \$\$\mathrm{logit}\\\lambda\_{t,d} =
\alpha_t + \beta d, \qquad d \in \\0, 1\\.\$\$ For each profile \\t\\
with at least one verified subject, \\\alpha_t\\ is the log odds of
verification for healthy subjects with profile \\t\\, and \\\alpha_t +
\beta\\ is the log odds of verification for diseased subjects with the
same profile. Profiles with no verified subject have \\\lambda\_{t,d} =
0\\ and no intercept.

\\\beta\\ is the log odds ratio of verification between diseased and
healthy subjects with the same profile: \$\$\beta = \log
\frac{\lambda\_{t,1}\\/\\(1 - \lambda\_{t,1})}{\lambda\_{t,0}\\/\\(1 -
\lambda\_{t,0})}.\$\$ Equivalently, \\\exp(-\beta)\\ is the odds ratio
of disease between unverified and verified subjects with the same
profile. When \\\beta = 0\\ verification does not depend on the disease
status given the profile, which is the missing at random assumption. A
positive \\\beta\\ means that diseased subjects are verified more often
than healthy subjects with the same profile, so that unverified subjects
are less often diseased than verified ones with the same profile; a
negative \\\beta\\ has the opposite meaning.

With \\\pi = P(D = 1)\\ and \\\theta\_{d,t} = P(T = t \mid D = d)\\ the
disease parameters, and \\n\_{t,d}\\, \\u_t\\, \\\pi_1 = \pi\\, \\\pi_0
= 1 - \pi\\ as in the counts of the observed data, the log-likelihood is
\$\$\ell = \sum_t \Big\[\sum_d n\_{t,d}
\log(\pi_d\\\theta\_{d,t}\\\lambda\_{t,d}) + u_t \log \sum_d
\pi_d\\\theta\_{d,t}\\(1 - \lambda\_{t,d})\Big\].\$\$

With a saturated disease model, which places no restriction on
\\\theta_d\\, the data identify the proportion of diseased among
verified subjects of each profile, but not \\\beta\\. The proportion of
diseased among unverified subjects with a verified profile \\t\\ is
\\\mathrm{expit}(\mathrm{logit}\\\phi_t - \beta)\\, where \\\phi_t\\ is
the proportion among verified subjects, so \\\beta\\ links the two
groups but the data say nothing about it: its posterior equals its
prior. The prior on \\\beta\\ is then an explicit assumption, and
varying it, or fixing \\\beta\\ at a range of values, shows how much the
results depend on it.

With the log-linear and random-effects disease models, \\\beta\\ is
identified through the restrictions that these models place on the
disease distribution. This information comes from the assumed structure,
not from the data alone, and it is fragile if the structure is wrong.

## See also

[`mar()`](https://giovannitinervia9.github.io/baysesp/reference/mar.md),
[verification_mechanism](https://giovannitinervia9.github.io/baysesp/reference/verification_mechanism.md),
[`sensitivity_analysis()`](https://giovannitinervia9.github.io/baysesp/reference/sensitivity_analysis.md)

## Examples

``` r
mnar()
#> <mnar> verification depends on the test results and on the disease
#>   odds_ratio_prior: Normal(mean = 0, sd = 1)
#>   intercepts_prior: Normal(mean = 0, sd = 2.5)
mnar(odds_ratio_prior = prior_normal(mean = log(2), sd = 0.5))
#> <mnar> verification depends on the test results and on the disease
#>   odds_ratio_prior: Normal(mean = 0.6931, sd = 0.5)
#>   intercepts_prior: Normal(mean = 0, sd = 2.5)
```
