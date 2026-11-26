# Accuracy of the individual tests

Posterior summaries of the prevalence and of the sensitivity (`Se`),
specificity (`Sp`), positive predictive value (`PPV`) and negative
predictive value (`NPV`) of each test (Begg and Greenes, 1983; Zhou,
1993; Zhou, 1994).

## Usage

``` r
test_accuracy(fit, level = 0.95)
```

## Arguments

- fit:

  A fitted model, from
  [`baysesp()`](https://giovannitinervia9.github.io/baysesp/reference/baysesp.md).

- level:

  Probability of the credible intervals.

## Value

A data frame with one row per quantity and columns `quantity`, `mean`,
`median`, `lower` and `upper` (equal-tailed credible interval) and
`identified`: `TRUE` if the quantity is estimated from the data alone,
`FALSE` if it also depends on the prior of a parameter that the data
cannot estimate (the proportion of diseased in never-verified profiles)
(Buzoianu and Kadane, 2008; Gustafson, 2015).

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

The first row is the prevalence \\\pi\\, the proportion of diseased
subjects in the population.

## References

Begg, C. B. and Greenes, R. A. (1983). Assessment of diagnostic tests
when disease verification is subject to selection bias. *Biometrics*,
39(1), 207–215. [doi:10.2307/2530820](https://doi.org/10.2307/2530820)

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

Zhou, X.-H. (1993). Maximum likelihood estimators of sensitivity and
specificity corrected for verification bias. *Communications in
Statistics – Theory and Methods*, 22(11), 3177–3198.
[doi:10.1080/03610929308831209](https://doi.org/10.1080/03610929308831209)

Zhou, X.-H. (1994). Effect of verification bias on positive and negative
predictive values. *Statistics in Medicine*, 13(17), 1737–1745.
[doi:10.1002/sim.4780131705](https://doi.org/10.1002/sim.4780131705)

## See also

[`rule_accuracy()`](https://giovannitinervia9.github.io/baysesp/reference/rule_accuracy.md)
for combinations of tests,
[`compare_tests()`](https://giovannitinervia9.github.io/baysesp/reference/compare_tests.md).

## Examples

``` r
# \donttest{
set.seed(1)
d <- simulate_study(n = 1000, prevalence = 0.2,
                    sensitivity = c(A = 0.9, B = 0.8),
                    specificity = c(A = 0.85, B = 0.95),
                    verification_prob = c(0.2, 0.9, 0.9, 1))
d <- baysesp_data(d, tests = c("A", "B"), verified = "verified", disease = "disease")
test_accuracy(baysesp(d, seed = 1), level = 0.9)
#>     quantity      mean    median     lower     upper identified
#> 1 prevalence 0.2032967 0.2029564 0.1818824 0.2254947       TRUE
#> 2      Se[A] 0.8949437 0.8973944 0.8506317 0.9317398       TRUE
#> 3      Sp[A] 0.8497024 0.8498352 0.8282673 0.8705420       TRUE
#> 4     PPV[A] 0.6030221 0.6031086 0.5549486 0.6504706       TRUE
#> 5     NPV[A] 0.9693508 0.9702629 0.9548735 0.9804861       TRUE
#> 6      Se[B] 0.7975647 0.7987924 0.7438641 0.8461167       TRUE
#> 7      Sp[B] 0.9472373 0.9474472 0.9336387 0.9597483       TRUE
#> 8     PPV[B] 0.7940863 0.7951280 0.7451814 0.8409058       TRUE
#> 9     NPV[B] 0.9482413 0.9488913 0.9328772 0.9616369       TRUE
# }
```
