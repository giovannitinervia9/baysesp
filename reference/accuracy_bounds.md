# Bounds for the absolute accuracy when some profiles are never verified

Gives the range of values that the prevalence and the sensitivity and
specificity of each test (or rule) can take, whatever the proportion of
diseased subjects in the response profiles that are never verified.

## Usage

``` r
accuracy_bounds(fit, rules = fit@data@tests, odds_ratio = NULL, level = 0.95)
```

## Arguments

- fit:

  A model fitted with
  [`baysesp()`](https://giovannitinervia9.github.io/baysesp/reference/baysesp.md)
  and
  [`saturated()`](https://giovannitinervia9.github.io/baysesp/reference/saturated.md).

- rules:

  Rules or test names, see
  [`rule_accuracy()`](https://giovannitinervia9.github.io/baysesp/reference/rule_accuracy.md).
  By default, the individual tests.

- odds_ratio:

  The odds ratio of verification between diseased and healthy subjects
  with the same response profile, \\\exp(\beta)\\: a single positive
  number (one value only; to compare several values, call the function
  once for each). 1 means that verification does not depend on the
  disease status (missing at random); values above 1 mean that the
  diseased are verified more often. The default, `NULL`, uses 1 for a
  model fitted with `verification = mar()`, and gives an error for a
  model fitted with `verification = mnar()`, for which a value must be
  given.

- level:

  Probability of the credible intervals.

## Value

A data frame with columns `quantity`, `lower_bound`, `upper_bound`
(posterior medians of the bounds), `lower` and `upper` (the lower
quantile of the lower bound and the upper quantile of the upper bound).

## Details

When some profiles (combinations of test results) are never verified,
the data say nothing about how many diseased subjects they contain, and
the absolute accuracy depends on it. A prior on that proportion gives a
single answer, which rests entirely on the prior (Buzoianu and Kadane,
2008). These bounds instead show what the data alone allow (Gustafson,
2015): the lowest and highest value of each quantity when the proportion
of diseased in the never-verified profiles ranges from 0 to 1. They do
not depend on the prior. Wide bounds mean that the absolute accuracy is
mostly determined by the assumption on the never-verified profiles.

The bounds also depend on the odds ratio of verification between
diseased and healthy subjects with the same profile, which the data
cannot estimate either: they are computed for one fixed value,
`odds_ratio` (by default 1, verification that does not depend on the
disease). Call the function again with other values to see how the
bounds move.

**Bounds**. Let \\\psi_t = P(D = 1 \mid T = t)\\ be the proportion of
diseased in profile \\t\\ and \\\omega_t = P(T = t)\\ its frequency. For
a profile with verified subjects \\\psi_t\\ is estimated from the data:
with odds ratio of verification \\OR = \exp(\beta)\\, \\\psi_t = v_t
\phi_t + (1 - v_t)\\\mathrm{expit}(\mathrm{logit}\\\phi_t - \beta)\\,
where \\v_t\\ is the fraction of verified subjects in the profile and
\\\phi_t\\ the proportion of diseased among them (with \\OR = 1\\,
\\\psi_t = \phi_t\\). For a never-verified profile \\\psi_t\\ is
unknown. For a rule \\R\\ that declares no never-verified profile
positive, let \$\$A_R = \sum\_{t: r_R(t) = 1} \omega_t \psi_t, \qquad
B_R = \sum\_{t: r_R(t) = 1} \omega_t (1 - \psi_t),\$\$ the proportions
of subjects who are declared positive and are diseased or healthy, \\c\\
the proportion of diseased subjects with a verified profile and
\\\omega_0\\ the frequency of the never-verified profiles. Then \\Se_R =
A_R / \pi\\ and \\1 - Sp_R = B_R / (1 - \pi)\\, and since \\c \le \pi
\le c + \omega_0\\, \$\$\frac{A_R}{c + \omega_0} \le Se_R \le
\frac{A_R}{c}, \qquad 1 - \frac{B_R}{1 - c - \omega_0} \le Sp_R \le 1 -
\frac{B_R}{1 - c}.\$\$ The upper bound of the sensitivity and the lower
bound of the specificity are the values of the naive analysis that
considers healthy every subject with a never-verified profile; the other
bounds correspond to considering all of them diseased.

The bounds are functions of identified parameters (and of the odds
ratio), so they have a posterior distribution: the result gives their
posterior medians and a credible region that contains the whole
identification interval with probability at least `level` (Gustafson,
2015). For rules that declare a never-verified profile positive the
bounds are `NA`.

The accuracy is defined for a *rule* \\R\\, which declares a subject
positive or negative from the response profile \\t\\, the vector of the
results of the tests: \\r_R(t) = 1\\ if the profile is declared
positive, 0 otherwise. A single test is the rule "this test is
positive". The quantities depend on the prevalence \\\pi = P(D = 1)\\
and on the distribution of the profiles in the diseased, \\\theta\_{1,t}
= P(T = t \mid D = 1)\\, and in the healthy, \\\theta\_{0,t} = P(T = t
\mid D = 0)\\:

- **sensitivity** (`Se`): \$\$Se_R = P(R \textrm{ positive} \mid D = 1)
  = \sum_t r_R(t)\\ \theta\_{1,t},\$\$ the proportion of diseased
  subjects that the rule declares positive. It is also called the true
  positive rate (TPR);

- **specificity** (`Sp`): \$\$Sp_R = P(R \textrm{ negative} \mid D = 0)
  = 1 - \sum_t r_R(t)\\ \theta\_{0,t},\$\$ the proportion of healthy
  subjects that the rule declares negative; \\1 - Sp_R\\ is the false
  positive rate (FPR);

- **positive predictive value** (`PPV`): \$\$PPV_R = P(D = 1 \mid R
  \textrm{ positive}) = \frac{\pi\\ Se_R}{\pi\\ Se_R + (1 - \pi)(1 -
  Sp_R)},\$\$ the probability that a subject declared positive is
  diseased;

- **negative predictive value** (`NPV`): \$\$NPV_R = P(D = 0 \mid R
  \textrm{ negative}) = \frac{(1 - \pi)\\ Sp_R}{(1 - \pi)\\ Sp_R + \pi\\
  (1 - Se_R)},\$\$ the probability that a subject declared negative is
  healthy.

Sensitivity and specificity describe the rule and do not depend on the
prevalence; the predictive values describe what a result means for a
subject, and change with the prevalence of the population. For combined
rules, such as
[`any_positive()`](https://giovannitinervia9.github.io/baysesp/reference/any_positive.md),
the accuracy depends on the joint distribution of the tests, not only on
their individual sensitivities and specificities (Torrance-Rynard and
Walter, 1997; Vacek, 1985).

Each quantity is computed from every posterior draw of \\\pi\\,
\\\theta_0\\ and \\\theta_1\\. The results report its posterior mean and
median and an equal-tailed credible interval with probability `level`.

## References

Buzoianu, M. and Kadane, J. B. (2008). Adjusting for verification bias
in diagnostic test evaluation: A Bayesian approach. *Statistics in
Medicine*, 27(13), 2453–2473.
[doi:10.1002/sim.3099](https://doi.org/10.1002/sim.3099)

Gustafson, P. (2015). *Bayesian Inference for Partially Identified
Models: Exploring the Limits of Limited Data*. Chapman & Hall/CRC, Boca
Raton, FL.

Torrance-Rynard, V. L. and Walter, S. D. (1997). Effects of dependent
errors in the assessment of diagnostic test performance. *Statistics in
Medicine*, 16(19), 2157–2175.
[doi:10.1002/(SICI)1097-0258(19971015)16:19\<2157::AID-SIM653\>3.0.CO;2-X](https://doi.org/10.1002/%28SICI%291097-0258%2819971015%2916%3A19%3C2157%3A%3AAID-SIM653%3E3.0.CO%3B2-X)

Vacek, P. M. (1985). The effect of conditional dependence on the
evaluation of diagnostic tests. *Biometrics*, 41(4), 959–968.
[doi:10.2307/2530967](https://doi.org/10.2307/2530967)

## See also

[`accuracy_ratios()`](https://giovannitinervia9.github.io/baysesp/reference/accuracy_ratios.md),
[`sensitivity_analysis()`](https://giovannitinervia9.github.io/baysesp/reference/sensitivity_analysis.md)

## Examples

``` r
# \donttest{
set.seed(1)
# profile (0,0) is never verified
d <- simulate_study(n = 2000, prevalence = 0.2,
                    sensitivity = c(A = 0.9, B = 0.8),
                    specificity = c(A = 0.85, B = 0.95),
                    verification_prob = c(0, 0.9, 0.9, 1))
d <- baysesp_data(d, tests = c("A", "B"), verified = "verified", disease = "disease")
fit <- baysesp(d, saturated(unverified_prevalence_prior = elicit_beta(0.03, 0.06)),
               seed = 1)
accuracy_bounds(fit)
#>     quantity lower_bound upper_bound     lower     upper
#> 1 prevalence   0.1998029   0.8397389 0.1825565 0.8555800
#> 2      Se[A]   0.2150756   0.9045453 0.1958355 0.9325346
#> 3      Sp[A]   0.2154108   0.8429136 0.1716957 0.8601760
#> 4      Se[B]   0.1912372   0.8044312 0.1726520 0.8420550
#> 5      Sp[B]   0.7360968   0.9472660 0.6848990 0.9580449
accuracy_bounds(fit, odds_ratio = 2)
#>     quantity lower_bound upper_bound     lower     upper
#> 1 prevalence   0.1965725   0.8365040 0.1795691 0.8525168
#> 2      Se[A]   0.2132777   0.9082051 0.1941214 0.9353522
#> 3      Sp[A]   0.2174917   0.8408455 0.1741221 0.8579767
#> 4      Se[B]   0.1906394   0.8123313 0.1721827 0.8492972
#> 5      Sp[B]   0.7348246   0.9461930 0.6850631 0.9568404
# }
```
