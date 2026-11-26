# Simulate a study with partial verification

Generates individual data for a study in which `K` binary tests are
applied to all subjects, but only some of them are verified by the
reference standard. For each subject the disease status is drawn first,
then the `K` test results and finally the verification, which depends on
the response profile (the vector of test results) and, when `odds_ratio`
is not 1, also on the disease status. Useful for examples, tests and
simulation studies.

## Usage

``` r
simulate_study(
  n,
  prevalence,
  sensitivity,
  specificity,
  verification_prob,
  odds_ratio = 1,
  dependence = independence()
)
```

## Arguments

- n:

  Number of subjects: a single whole number, at least 1. No default.

- prevalence:

  Probability of disease: a single number in \\\[0, 1\]\\. No default.

- sensitivity:

  Sensitivity of each test: a numeric vector with values in \\\[0,
  1\]\\, one element per test, of length `K >= 1`. No default. The names
  of `sensitivity`, if present, become the names of the test columns of
  the output; if it has no names, the tests are called `T1`, ..., `TK`.
  The names `verified` and `disease` are not allowed.

- specificity:

  Specificity of each test: a numeric vector with values in \\\[0, 1\]\\
  and the same length as `sensitivity`; its k-th element is the
  specificity of the k-th test. Its names are ignored. No default.

- verification_prob:

  Probability of verification of healthy subjects, with values in \\\[0,
  1\]\\. It is either a single number, which applies to all response
  profiles, or a vector of length \\2^K\\ with one number per response
  profile, in binary order: the profile \\t\\ is at position \\1 +
  \sum_k t_k 2^{k-1}\\, so the first profile is the one with all tests
  negative and the first test changes fastest. With `K = 2` the order is
  (0,0), (1,0), (0,1), (1,1). A value 0 means that the profile is never
  verified and a value 1 that it is always verified. No default.

- odds_ratio:

  Odds ratio of verification between diseased and healthy subjects with
  the same response profile: a single positive number, the same for all
  profiles. Default `1`, which gives verification that is missing at
  random (it does not depend on the disease status). Values greater than
  1 make diseased subjects more likely to be verified, values between 0
  and 1 less likely; any value different from 1 gives verification that
  is missing not at random.

- dependence:

  How the test results depend on each other given the disease status: a
  single dependence mechanism, one of `dependence = independence()` (the
  tests are independent given the disease status),
  `dependence = latent_dependence(...)` (the tests share a normal latent
  variable, with one loading per test and disease class) or
  `dependence = probit_dependence(...)` (the tests have jointly normal
  latent variables, with one correlation per pair of tests and disease
  class). Default
  [`independence()`](https://giovannitinervia9.github.io/baysesp/reference/independence.md).
  Whatever the mechanism, `sensitivity` and `specificity` remain the
  marginal sensitivity and specificity of each test.

## Value

A data frame with one row per subject and the following columns: one 0/1
column for each test, always observed, named as described in
`sensitivity`; `verified` (0/1); `disease` (0/1 for verified subjects,
`NA` for unverified subjects).

## Details

For each subject the function simulates:

- the disease status \\D\\, Bernoulli with probability `prevalence`;

- the `K` test results \\T_1, \dots, T_K\\ given \\D\\, with \\P(T_k = 1
  \mid D = 1)\\ equal to `sensitivity[k]` and \\P(T_k = 0 \mid D = 0)\\
  equal to `specificity[k]`. The tests are independent given \\D\\ with
  `dependence = independence()`, and dependent given \\D\\ with the
  other values of `dependence` (see the paragraph on the dependence
  between tests below);

- the verification \\V\\, Bernoulli with probability \\\lambda\_{t,d} =
  P(V = 1 \mid T = t, D = d)\\.

The verification probability of a healthy subject with profile \\t\\ is
\\\lambda\_{t,0}\\ = `verification_prob` for that profile. The
verification probability of a diseased subject with the same profile is
\\\lambda\_{t,1} = \mathrm{expit}(\mathrm{logit}\\\lambda\_{t,0} +
\log(\mathrm{odds\\ratio}))\\, so that the odds of verification of
diseased subjects are `odds_ratio` times the odds for healthy subjects
with the same profile.

**Dependence between the tests.** The argument `dependence` sets the
joint distribution of the `K` results within each disease class, that is
the probabilities \\\theta\_{d,t} = P(T = t \mid D = d)\\ of the
response profiles \\t\\ among the diseased (\\d = 1\\) and among the
healthy (\\d = 0\\). Every mechanism keeps the sensitivities and
specificities given in `sensitivity` and `specificity` as the marginal
probabilities, and changes only how often the tests are right or wrong
together:

- `dependence = independence()`: \\\theta\_{1,t} = \prod_k
  Se_k^{t_k}(1 - Se_k)^{1 - t_k}\\ and \\\theta\_{0,t} = \prod_k (1 -
  Sp_k)^{t_k} Sp_k^{1 - t_k}\\;

- `dependence = latent_dependence(diseased_loadings, healthy_loadings)`:
  each subject has a latent variable \\Z \sim N(0, 1)\\, independent of
  \\D\\, and given \\D = d\\ and \\Z = z\\ the tests are independent
  with \\P(T_k = 1 \mid D = d, Z = z) = \Phi(a\_{dk} + b\_{dk} z)\\,
  where \\b\_{dk}\\ is the loading of test \\k\\ in class \\d\\ and
  \\a\_{dk}\\ is chosen so that the marginal probabilities are the
  sensitivities and specificities. Tests \\j\\ and \\k\\ are positively
  dependent when \\b\_{dj} b\_{dk} \> 0\\ and negatively dependent when
  \\b\_{dj} b\_{dk} \< 0\\;

- `dependence = probit_dependence(diseased_correlation, healthy_correlation)`:
  given \\D = d\\, test \\k\\ is positive when a standard normal latent
  variable \\Y_k\\ falls below a threshold that gives the required
  sensitivity or specificity, and the latent variables of the `K` tests
  are jointly normal with correlation matrix \\R_d\\. The element
  \\R\_{d,jk}\\ is the correlation, of any sign, between tests \\j\\ and
  \\k\\ in class \\d\\.

The size of the loadings or of the correlation matrices is checked here
against the number of tests `K`, which the dependence mechanism does not
know when it is created. A correlation matrix that is not positive
definite for `K` tests gives an error.

**Supported data-generating mechanisms.**

- any number `K >= 1` of tests, with arbitrary prevalence, sensitivities
  and specificities, with tests that are independent given the disease
  status or dependent through one of the mechanisms above;

- verification that depends on the response profile only, with the same
  probability for all profiles (a single `verification_prob`) or with
  one probability per profile. A value 0 means that the profile is never
  verified, a value 1 that it is always verified (for healthy subjects;
  see `odds_ratio` for diseased subjects);

- missing at random verification (`odds_ratio = 1`), where verification
  depends only on the profile, and missing not at random verification
  (`odds_ratio` different from 1), where the odds of verification also
  depend on the disease status through the same odds ratio for all
  profiles. Profiles with verification probability 0 or 1 are never
  verified, or always verified, whatever `odds_ratio` is.

**Not supported:** dependence between the tests given the disease status
other than through the two latent mechanisms above (for example
log-linear interactions given as odds ratios), covariates, verification
that depends on something other than the response profile and the
disease status, and a disease-dependent verification whose odds ratio
differs across profiles.

## See also

[`independence()`](https://giovannitinervia9.github.io/baysesp/reference/independence.md),
[`latent_dependence()`](https://giovannitinervia9.github.io/baysesp/reference/latent_dependence.md),
[`probit_dependence()`](https://giovannitinervia9.github.io/baysesp/reference/probit_dependence.md)

## Examples

``` r
set.seed(1)

# Two tests, verification depending only on the results (MAR): subjects
# negative to both tests are verified with probability 0.1, the others with
# probability 0.9 or 1. The profiles are in the order (0,0), (1,0), (0,1), (1,1).
d <- simulate_study(n = 500, prevalence = 0.2,
                    sensitivity = c(A = 0.9, B = 0.7),
                    specificity = c(A = 0.85, B = 0.95),
                    verification_prob = c(0.1, 0.9, 0.9, 1))
head(d)
#>   A B verified disease
#> 1 0 0        0      NA
#> 2 0 0        0      NA
#> 3 0 0        0      NA
#> 4 1 0        1       1
#> 5 0 0        0      NA
#> 6 0 1        1       1
# proportion of verified subjects by profile
tapply(d$verified, paste(d$A, d$B), mean)
#>       0 0       0 1       1 0       1 1 
#> 0.1057402 0.8333333 0.9078947 1.0000000 

# Everyone verified: a single verification probability equal to 1
d_full <- simulate_study(n = 500, prevalence = 0.2,
                         sensitivity = c(A = 0.9, B = 0.7),
                         specificity = c(A = 0.85, B = 0.95),
                         verification_prob = 1)
mean(d_full$verified)
#> [1] 1

# A never-verified profile: only subjects positive to at least one test are
# verified, so the disease status is missing for every (0,0) subject
d_screen <- simulate_study(n = 500, prevalence = 0.2,
                           sensitivity = c(A = 0.9, B = 0.7),
                           specificity = c(A = 0.85, B = 0.95),
                           verification_prob = c(0, 1, 1, 1))
table(paste(d_screen$A, d_screen$B), d_screen$verified)
#>      
#>         0   1
#>   0 0 305   0
#>   0 1   0  23
#>   1 0   0  99
#>   1 1   0  73

# Verification that depends on the disease (MNAR): with odds_ratio = 3 the
# odds of verification are three times higher for the diseased, so among
# the verified the diseased are over-represented
d_mnar <- simulate_study(n = 5000, prevalence = 0.2,
                         sensitivity = c(A = 0.9, B = 0.7),
                         specificity = c(A = 0.85, B = 0.95),
                         verification_prob = 0.3, odds_ratio = 3)
mean(d_mnar$disease, na.rm = TRUE)   # well above the prevalence 0.2
#> [1] 0.3348066

# Three tests without names: the columns are called T1, T2, T3, and there
# are 2^3 = 8 profiles, so verification_prob has 8 values
d3 <- simulate_study(n = 500, prevalence = 0.1,
                     sensitivity = c(0.9, 0.8, 0.7),
                     specificity = c(0.9, 0.9, 0.95),
                     verification_prob = c(0.05, rep(0.8, 6), 1))
names(d3)
#> [1] "T1"       "T2"       "T3"       "verified" "disease" 

# Three tests dependent given the disease: among the diseased tests A and B
# share a latent severity, C much less; among the healthy only A and B
# are (weakly) dependent
d_dep <- simulate_study(n = 5000, prevalence = 0.2,
                        sensitivity = c(A = 0.9, B = 0.8, C = 0.7),
                        specificity = c(A = 0.85, B = 0.9, C = 0.95),
                        verification_prob = c(0.05, rep(0.8, 6), 1),
                        dependence = latent_dependence(
                          diseased_loadings = c(1, 1, 0.3),
                          healthy_loadings = c(0.5, 0.5, 0)))
sick <- which(d_dep$disease == 1)
cor(d_dep[sick, c("A", "B", "C")])
#>            A          B          C
#> A 1.00000000 0.10113242 0.02708195
#> B 0.10113242 1.00000000 0.07264363
#> C 0.02708195 0.07264363 1.00000000

# Correlation 0.4 between every pair of tests among the diseased and -0.2
# among the healthy
d_probit <- simulate_study(n = 5000, prevalence = 0.2,
                           sensitivity = c(A = 0.9, B = 0.8, C = 0.7),
                           specificity = c(A = 0.85, B = 0.9, C = 0.95),
                           verification_prob = 1,
                           dependence = probit_dependence(
                             diseased_correlation = 0.4,
                             healthy_correlation = -0.2))
```
