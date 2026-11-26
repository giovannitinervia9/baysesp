# Ratios between the accuracy of tests when some profiles are never verified

Gives, for each pair of tests (or rules), the ratio of their
sensitivities and the ratio of their false positive rates (one minus the
specificity). Unlike the absolute accuracy, these ratios are estimated
by the data even when some response profiles are never verified (Alonzo,
2005; Pepe and Alonzo, 2001; Schatzkin et al., 1987).

## Usage

``` r
accuracy_ratios(fit, rules = fit@data@tests, odds_ratio = NULL, level = 0.95)
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

A data frame with one row per ratio and pair of rules and columns
`quantity`, `mean`, `median`, `lower`, `upper` and `probability`, the
posterior probability that the ratio is greater than 1: for rTPR, that
the first rule is more sensitive; for rFPR, that it is less specific.
Only the rules that declare no never-verified profile positive are
compared.

## Details

When some profiles (combinations of test results) are never verified,
the data say nothing about how many diseased subjects they contain
(Gustafson, 2015), and the sensitivities and specificities depend on it.
If neither test declares a never-verified profile positive, however, the
unknown diseased subjects are missed by both tests: they lower the two
sensitivities by the same factor, and the ratio of the sensitivities
does not change. The same holds for the false positive rates. These
ratios therefore tell how the tests compare, for example that one
detects 15% more diseased subjects than the other, without any
assumption on the never-verified profiles; the posterior probability
that a ratio is greater than 1 tells which test is more sensitive, or
less specific. Tests or rules that declare a never-verified profile
positive are left out.

The ratios depend on the odds ratio of verification between diseased and
healthy subjects with the same profile, which the data cannot estimate:
they are computed for one fixed value, `odds_ratio` (by default 1,
verification that does not depend on the disease).

**Why the ratios are estimated by the data.** Let \\\omega_t\\ be the
frequency of profile \\t\\, \\\psi_t\\ the proportion of diseased in it,
known from the data for the profiles with verified subjects, and \\\pi\\
the prevalence, which depends on the unknown \\\psi_t\\ of the
never-verified profiles. For a rule \\R\\ that declares no
never-verified profile positive, \$\$A_R = \sum\_{t: r_R(t) = 1}
\omega_t \psi_t, \qquad B_R = \sum\_{t: r_R(t) = 1} \omega_t (1 -
\psi_t)\$\$ are known, and \\Se_R = A_R / \pi\\, \\1 - Sp_R = B_R / (1 -
\pi)\\. For two such rules the unknown prevalence cancels out in the
ratios \$\$rTPR = \frac{Se_R}{Se\_{R'}} = \frac{A_R}{A\_{R'}}, \qquad
rFPR = \frac{1 - Sp_R}{1 - Sp\_{R'}} = \frac{B_R}{B\_{R'}},\$\$ which
are therefore estimated from the data alone (Alonzo, 2005; Pepe and
Alonzo, 2001; Schatzkin et al., 1987). \\rTPR \> 1\\ means that \\R\\ is
more sensitive than \\R'\\, \\rFPR \> 1\\ that it is less specific,
whatever the proportion of diseased in the never-verified profiles.

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

Alonzo, T. A. (2005). Verification bias-corrected estimators of the
relative true and false positive rates of two binary screening tests.
*Statistics in Medicine*, 24(3), 403–417.
[doi:10.1002/sim.1959](https://doi.org/10.1002/sim.1959)

Gustafson, P. (2015). *Bayesian Inference for Partially Identified
Models: Exploring the Limits of Limited Data*. Chapman & Hall/CRC, Boca
Raton, FL.

Pepe, M. S. and Alonzo, T. A. (2001). Comparing disease screening tests
when true disease status is ascertained only for screen positives.
*Biostatistics*, 2(3), 249–260.
[doi:10.1093/biostatistics/2.3.249](https://doi.org/10.1093/biostatistics/2.3.249)

Schatzkin, A., Connor, R. J., Taylor, P. R. and Bunnag, B. (1987).
Comparing new and old screening tests when a reference procedure cannot
be performed on all screenees. *American Journal of Epidemiology*,
125(4), 672–678.
[doi:10.1093/oxfordjournals.aje.a114580](https://doi.org/10.1093/oxfordjournals.aje.a114580)

Torrance-Rynard, V. L. and Walter, S. D. (1997). Effects of dependent
errors in the assessment of diagnostic test performance. *Statistics in
Medicine*, 16(19), 2157–2175.
[doi:10.1002/(SICI)1097-0258(19971015)16:19\<2157::AID-SIM653\>3.0.CO;2-X](https://doi.org/10.1002/%28SICI%291097-0258%2819971015%2916%3A19%3C2157%3A%3AAID-SIM653%3E3.0.CO%3B2-X)

Vacek, P. M. (1985). The effect of conditional dependence on the
evaluation of diagnostic tests. *Biometrics*, 41(4), 959–968.
[doi:10.2307/2530967](https://doi.org/10.2307/2530967)

## See also

[`accuracy_bounds()`](https://giovannitinervia9.github.io/baysesp/reference/accuracy_bounds.md),
[`compare_tests()`](https://giovannitinervia9.github.io/baysesp/reference/compare_tests.md)

## Examples

``` r
# \donttest{
set.seed(1)
d <- simulate_study(n = 2000, prevalence = 0.2,
                    sensitivity = c(A = 0.9, B = 0.8),
                    specificity = c(A = 0.85, B = 0.95),
                    verification_prob = c(0, 0.9, 0.9, 1))
d <- baysesp_data(d, tests = c("A", "B"), verified = "verified", disease = "disease")
fit <- baysesp(d, saturated(unverified_prevalence_prior = elicit_beta(0.03, 0.06)),
               seed = 1)
accuracy_ratios(fit, rules = list("A", "B", any_positive()))
#>                     quantity      mean    median     lower     upper
#> 1            TPR[A] / TPR[B] 1.1257685 1.1243562 1.0559413 1.2038481
#> 2            FPR[A] / FPR[B] 2.9924388 2.9679967 2.3650247 3.7902284
#> 3 TPR[A] / TPR[any positive] 0.9037443 0.9044398 0.8722572 0.9318102
#> 4 FPR[A] / FPR[any positive] 0.7846215 0.7847116 0.7381380 0.8285368
#> 5 TPR[B] / TPR[any positive] 0.8033899 0.8041101 0.7607886 0.8432921
#> 6 FPR[B] / FPR[any positive] 0.2651496 0.2645190 0.2176734 0.3149205
#>   probability
#> 1           1
#> 2           1
#> 3           0
#> 4           0
#> 5           0
#> 6           0
# }
```
