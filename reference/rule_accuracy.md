# Accuracy of diagnostic rules

Posterior summaries of the sensitivity (`Se`), specificity (`Sp`),
positive predictive value (`PPV`) and negative predictive value (`NPV`)
of rules that combine the tests (Begg and Greenes, 1983; Zhou, 1993;
Zhou, 1994). For a single test the result is the same as in
[`test_accuracy()`](https://giovannitinervia9.github.io/baysesp/reference/test_accuracy.md).

## Usage

``` r
rule_accuracy(fit, rules, level = 0.95)
```

## Arguments

- fit:

  A fitted model, from
  [`baysesp()`](https://giovannitinervia9.github.io/baysesp/reference/baysesp.md).

- rules:

  A rule, such as
  [`any_positive()`](https://giovannitinervia9.github.io/baysesp/reference/any_positive.md),
  the name of a test, or a list of them.

- level:

  Probability of the credible intervals.

## Value

A data frame as in
[`test_accuracy()`](https://giovannitinervia9.github.io/baysesp/reference/test_accuracy.md),
without the row of the prevalence.

## Details

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

Begg, C. B. and Greenes, R. A. (1983). Assessment of diagnostic tests
when disease verification is subject to selection bias. *Biometrics*,
39(1), 207–215. [doi:10.2307/2530820](https://doi.org/10.2307/2530820)

Torrance-Rynard, V. L. and Walter, S. D. (1997). Effects of dependent
errors in the assessment of diagnostic test performance. *Statistics in
Medicine*, 16(19), 2157–2175.
[doi:10.1002/(SICI)1097-0258(19971015)16:19\<2157::AID-SIM653\>3.0.CO;2-X](https://doi.org/10.1002/%28SICI%291097-0258%2819971015%2916%3A19%3C2157%3A%3AAID-SIM653%3E3.0.CO%3B2-X)

Vacek, P. M. (1985). The effect of conditional dependence on the
evaluation of diagnostic tests. *Biometrics*, 41(4), 959–968.
[doi:10.2307/2530967](https://doi.org/10.2307/2530967)

Zhou, X.-H. (1993). Maximum likelihood estimators of sensitivity and
specificity corrected for verification bias. *Communications in
Statistics – Theory and Methods*, 22(11), 3177–3198.
[doi:10.1080/03610929308831209](https://doi.org/10.1080/03610929308831209)

Zhou, X.-H. (1994). Effect of verification bias on positive and negative
predictive values. *Statistics in Medicine*, 13(17), 1737–1745.
[doi:10.1002/sim.4780131705](https://doi.org/10.1002/sim.4780131705)

## Examples

``` r
# \donttest{
set.seed(1)
d <- simulate_study(n = 1000, prevalence = 0.2,
                    sensitivity = c(A = 0.9, B = 0.8),
                    specificity = c(A = 0.85, B = 0.95),
                    verification_prob = c(0.2, 0.9, 0.9, 1))
d <- baysesp_data(d, tests = c("A", "B"), verified = "verified", disease = "disease")
fit <- baysesp(d, seed = 1)
rule_accuracy(fit, list(any_positive(), all_positive()))
#>            quantity      mean    median     lower     upper identified
#> 1  Se[any positive] 0.9891497 0.9953785 0.9468068 0.9999939       TRUE
#> 2  Sp[any positive] 0.8037203 0.8040117 0.7748167 0.8314678       TRUE
#> 3 PPV[any positive] 0.5624813 0.5624924 0.5081414 0.6151883       TRUE
#> 4 NPV[any positive] 0.9965089 0.9985393 0.9822992 0.9999981       TRUE
#> 5  Se[all positive] 0.7033587 0.7042886 0.6338550 0.7687375       TRUE
#> 6  Sp[all positive] 0.9932194 0.9936070 0.9867010 0.9975990       TRUE
#> 7 PPV[all positive] 0.9635747 0.9655511 0.9297147 0.9869545       TRUE
#> 8 NPV[all positive] 0.9291343 0.9296021 0.9086617 0.9464525       TRUE
# }
```
