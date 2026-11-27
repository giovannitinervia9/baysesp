# Missing not at random verification

`mnar()` is the assumption that verification depends on the test results
and also on the disease status (Kosinski and Barnhart, 2003a): among
subjects with the same profile, diseased and healthy subjects can be
verified with different probabilities. Its arguments are the prior on
the log odds ratio of verification, the prior on the verification
intercepts and the choice of whether the odds ratio varies across
profiles.

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
  is not accepted. The default `prior_normal(0, 1)` centers \\\beta\\ at
  0, which is missing at random, and gives 95% of its prior mass to odds
  ratios between about 0.14 and 7.1.

- intercepts_prior:

  Prior on the verification intercepts \\\alpha_t\\, one for each
  profile with at least one verified subject. It must act on the real
  line, as for `odds_ratio_prior`. It is used only by the log-linear and
  the random-effects models; the saturated model does not use it. The
  default is `prior_normal(0, 2.5)`. In profiles where almost all
  subjects are diseased, the data say little about the verification
  probability of the healthy, \\\mathrm{expit}(\alpha_t)\\, and the
  prior pulls it towards 0.5; the shared \\\beta\\ then absorbs the
  difference. As a result the estimated odds ratio depends on this
  prior: in simulations with verification probabilities of 0.8,
  `prior_normal(0, 1.5)` gave a posterior median of the odds ratio of
  about 2.4 when the true value was 1, the default about 1.7. The
  sensitivities and specificities changed by less than 0.005. Read the
  odds ratio of the structured models together with
  [`sensitivity_analysis()`](https://giovannitinervia9.github.io/baysesp/reference/sensitivity_analysis.md).

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
0\\ and no intercept; in the log-linear and random effects models
[`baysesp()`](https://giovannitinervia9.github.io/baysesp/reference/baysesp.md)
warns when such a profile has subjects, because for it verification is
then not missing not at random.

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
groups but the data say nothing about it: its posterior equals its prior
(Gustafson, 2015). The prior on \\\beta\\ is then an explicit
assumption, and varying it, or fixing \\\beta\\ at a range of values,
shows how much the results depend on it (Buzoianu and Kadane, 2008;
Kosinski and Barnhart, 2003b).

With the log-linear (Baker, 1995) and random-effects (Dendukuri and
Joseph, 2001; Qu et al., 1996) disease models, \\\beta\\ is identified
through the restrictions that these models place on the disease
distribution. This information comes from the assumed structure, not
from the data alone, and it is fragile if the structure is wrong (Albert
and Dodd, 2004).

## References

Albert, P. S. and Dodd, L. E. (2004). A cautionary note on the
robustness of latent class models for estimating diagnostic error
without a gold standard. *Biometrics*, 60(2), 427–435.
[doi:10.1111/j.0006-341X.2004.00187.x](https://doi.org/10.1111/j.0006-341X.2004.00187.x)

Baker, S. G. (1995). Evaluating multiple diagnostic tests with partial
verification. *Biometrics*, 51(1), 330–337.
[doi:10.2307/2533339](https://doi.org/10.2307/2533339)

Buzoianu, M. and Kadane, J. B. (2008). Adjusting for verification bias
in diagnostic test evaluation: A Bayesian approach. *Statistics in
Medicine*, 27(13), 2453–2473.
[doi:10.1002/sim.3099](https://doi.org/10.1002/sim.3099)

Dendukuri, N. and Joseph, L. (2001). Bayesian approaches to modeling the
conditional dependence between multiple diagnostic tests. *Biometrics*,
57(1), 158–167.
[doi:10.1111/j.0006-341X.2001.00158.x](https://doi.org/10.1111/j.0006-341X.2001.00158.x)

Gustafson, P. (2015). *Bayesian Inference for Partially Identified
Models: Exploring the Limits of Limited Data*. Chapman & Hall/CRC, Boca
Raton, FL.

Kosinski, A. S. and Barnhart, H. X. (2003a). Accounting for nonignorable
verification bias in assessment of diagnostic tests. *Biometrics*,
59(1), 163–171.
[doi:10.1111/1541-0420.00019](https://doi.org/10.1111/1541-0420.00019)

Kosinski, A. S. and Barnhart, H. X. (2003b). A global sensitivity
analysis of performance of a medical diagnostic test when verification
bias is present. *Statistics in Medicine*, 22(17), 2711–2721.
[doi:10.1002/sim.1517](https://doi.org/10.1002/sim.1517)

Qu, Y., Tan, M. and Kutner, M. H. (1996). Random effects models in
latent class analysis for evaluating accuracy of diagnostic tests.
*Biometrics*, 52(3), 797–810.
[doi:10.2307/2533043](https://doi.org/10.2307/2533043)

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
