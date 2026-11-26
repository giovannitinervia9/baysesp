# Saturated model

The model with no assumption on the dependence between the tests: the
distribution of the profiles is free within each disease class. It is
written in a transparent parametrization, which separates the identified
parameters from the non-identified ones:

- identified: the frequency of each response profile, the verification
  probability of each profile with at least one verified subject, and
  the proportion of diseased among its verified subjects;

- not identified: the proportion of diseased among the unverified
  subjects of each such profile, which is fixed by the verification
  mechanism, and the proportion of diseased in each never-verified
  profile, whose posterior equals its prior.

## Usage

``` r
saturated(
  frequency_prior = prior_dirichlet(1),
  verified_prior = prior_beta(0.5, 2),
  profile_prevalence_prior = prior_beta(0.5, 2),
  unverified_prevalence_prior = NULL,
  verification = mar()
)
```

## Arguments

- frequency_prior:

  [`prior_dirichlet()`](https://giovannitinervia9.github.io/baysesp/reference/prior_dirichlet.md)
  on the vector of profile frequencies \\\omega = (\omega_1, \dots,
  \omega_J)\\. Its `concentration` is a single positive number, the same
  for every profile, or a vector with one positive value per profile, in
  binary order: the profile with all tests negative comes first, and the
  first test changes fastest (with two tests: (0,0), (1,0), (0,1),
  (1,1)). Default `prior_dirichlet(1)`, the uniform distribution on the
  simplex.

- verified_prior:

  Prior on the verification probability \\v_t\\ of each profile with at
  least one verified subject. The same prior is used, independently, for
  each such profile. A prior on a probability, such as
  [`prior_beta()`](https://giovannitinervia9.github.io/baysesp/reference/prior_beta.md)
  or
  [`elicit_beta()`](https://giovannitinervia9.github.io/baysesp/reference/elicit_beta.md),
  acts on \\v_t\\ itself; a prior on the real line, such as
  [`prior_normal()`](https://giovannitinervia9.github.io/baysesp/reference/prior_normal.md)
  or
  [`prior_student()`](https://giovannitinervia9.github.io/baysesp/reference/prior_student.md),
  acts on \\\mathrm{logit}\\v_t\\. Default `prior_beta(0.5, 2)`, the
  uniform distribution on (0, 1). Never-verified profiles have no such
  parameter.

- profile_prevalence_prior:

  Prior on \\\phi^{(1)}\_t\\, the proportion of diseased among the
  verified subjects of each profile with at least one verified subject.
  The same prior is used, independently, for each such profile, with the
  same rules as `verified_prior`: a prior on a probability acts on the
  proportion, a prior on the real line on its logit. Default
  `prior_beta(0.5, 2)`, the uniform distribution on (0, 1).

- unverified_prevalence_prior:

  Prior on \\q_t\\, the proportion of diseased in each never-verified
  profile (a profile with no verified subject). It has no default prior:
  its default `NULL` means that no prior is given. It is required when
  at least one profile is never verified, and
  [`baysesp()`](https://giovannitinervia9.github.io/baysesp/reference/baysesp.md)
  stops with an error if it is missing in that case. Its value can be
  one prior, which applies to every never-verified profile, or a list
  with one prior per never-verified profile, in binary order among those
  profiles (with two tests, if both (0,1) and (1,1) are never verified,
  the list has two elements: the first for (0,1), the second for (1,1)).
  A list of length one is the same as a single prior. Each prior has the
  same allowed classes as `verified_prior`. Its posterior equals its
  prior, so it is an explicit assumption; use
  [`elicit_beta()`](https://giovannitinervia9.github.io/baysesp/reference/elicit_beta.md)
  to specify it from a plausible value and an upper bound.

- verification:

  Verification mechanism:
  [`mar()`](https://giovannitinervia9.github.io/baysesp/reference/mar.md)
  (the default), where verification depends only on the test results, or
  [`mnar()`](https://giovannitinervia9.github.io/baysesp/reference/mnar.md),
  where it also depends on the disease status. See Details for what
  changes.

## Value

An object of class `saturated`.

## Details

With \\K\\ tests there are \\J = 2^K\\ response profiles \\t\\. For each
profile the model has four parameters:

- \\\omega_t = P(T = t)\\: frequency of the profile among all subjects;

- \\v_t = P(V = 1 \mid T = t)\\: probability that a subject with profile
  \\t\\ is verified;

- \\\phi^{(1)}\_t = P(D = 1 \mid T = t, V = 1)\\: proportion of diseased
  among the verified subjects with profile \\t\\;

- \\\phi^{(0)}\_t = P(D = 1 \mid T = t, V = 0)\\: proportion of diseased
  among the unverified subjects with profile \\t\\.

Let \\n\_{t,1}\\ and \\n\_{t,0}\\ be the numbers of verified diseased
and verified healthy subjects with profile \\t\\, \\u_t\\ the number of
unverified subjects, and \\N_t = n\_{t,1} + n\_{t,0} + u_t\\. The
log-likelihood is \$\$\ell = \sum_t \Big\[N_t \log \omega_t +
(n\_{t,1} + n\_{t,0}) \log v_t + u_t \log(1 - v_t) + n\_{t,1} \log
\phi^{(1)}\_t + n\_{t,0} \log(1 - \phi^{(1)}\_t)\Big\].\$\$ The
unverified subjects enter only through \\N_t\\ and \\u_t\\, because
their disease status is not observed. The proportion \\\phi^{(0)}\_t\\
does not enter the log-likelihood: the data carry no information on it,
so it is determined by the verification mechanism and by
\\\phi^{(1)}\_t\\. A never-verified profile has no verified subject, so
\\v_t = 0\\ and only its first term remains.

With `verification = mar()` (missing at random) verification depends
only on the test results: among subjects with the same profile, verified
and unverified subjects have the same proportion of diseased, so
\\\phi^{(0)}\_t = \phi^{(1)}\_t\\.

With `verification = mnar()` (missing not at random) verification may
also depend on the disease status: \$\$\mathrm{logit}\\\phi^{(0)}\_t =
\mathrm{logit}\\\phi^{(1)}\_t - \beta,\$\$ where \\\beta\\ is the log
odds ratio of verification between diseased and healthy subjects with
the same profile. Equivalently, \\\exp(-\beta)\\ is the odds ratio of
disease between unverified and verified subjects with the same profile:
\\\beta \> 0\\ means that the unverified are less often diseased. There
is one \\\beta\\ common to all profiles with verified subjects, or one
per profile with `verification = mnar(by_profile = TRUE)`. Since
\\\beta\\ does not enter the log-likelihood, its posterior equals its
prior, and the results depend on the assumed value: they should be
examined over a range of plausible values.

The quantities of interest are obtained as follows. For each profile,
\$\$\psi_t = v_t \phi^{(1)}\_t + (1 - v_t) \phi^{(0)}\_t\$\$ is the
proportion of diseased among subjects with profile \\t\\; for a
never-verified profile \\\phi^{(0)}\_t = q_t\\ is the proportion of
diseased in it, so \\\psi_t = q_t\\. Then \$\$\pi = \sum_t \omega_t
\psi_t, \qquad \theta\_{1,t} = \frac{\omega_t \psi_t}{\pi}, \qquad
\theta\_{0,t} = \frac{\omega_t (1 - \psi_t)}{1 - \pi},\$\$ where \\\pi\\
is the prevalence, \\\theta\_{1,t} = P(T = t \mid D = 1)\\ and
\\\theta\_{0,t} = P(T = t \mid D = 0)\\. The sensitivity of a test is
the sum of \\\theta\_{1,t}\\ over the profiles with a positive result
for that test, and its specificity the sum of \\\theta\_{0,t}\\ over the
profiles with a negative result. Under both mechanisms the
never-verified profiles enter \\\pi\\ through \\q_t\\, whose posterior
equals its prior.

## Examples

``` r
saturated()
#> <saturated> saturated model, transparent parametrization
#>   frequency_prior:             Dirichlet(concentration = 1)
#>   verified_prior:              Beta(mean = 0.5, precision = 2)
#>   profile_prevalence_prior:    Beta(mean = 0.5, precision = 2)
#>   unverified_prevalence_prior: none
#>   verification:                mar
saturated(unverified_prevalence_prior = elicit_beta(mean = 0.005, upper = 0.015))
#> <saturated> saturated model, transparent parametrization
#>   frequency_prior:             Dirichlet(concentration = 1)
#>   verified_prior:              Beta(mean = 0.5, precision = 2)
#>   profile_prevalence_prior:    Beta(mean = 0.5, precision = 2)
#>   unverified_prevalence_prior: Beta(mean = 0.005, precision = 197.8)
#>   verification:                mar
saturated(verification = mnar())
#> <saturated> saturated model, transparent parametrization
#>   frequency_prior:             Dirichlet(concentration = 1)
#>   verified_prior:              Beta(mean = 0.5, precision = 2)
#>   profile_prevalence_prior:    Beta(mean = 0.5, precision = 2)
#>   unverified_prevalence_prior: none
#>   verification:                mnar
```
