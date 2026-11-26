# Sensitivity analysis for the assumptions the data cannot check

Shows how the accuracy of the tests changes with the two assumptions of
the saturated model that the data cannot check: the odds ratio of
verification between diseased and healthy subjects (set with
`verification = mnar()` in the model) and the prior on the proportion of
diseased in the never-verified profiles (Buzoianu and Kadane, 2008). For
the saturated model no new fit is needed: these parameters do not enter
the likelihood (Gustafson, 2015), so the draws of the identified
parameters are reused and only the assumptions are changed.

## Usage

``` r
sensitivity_analysis(
  fit,
  odds_ratios = NULL,
  unverified_priors = NULL,
  rules = NULL,
  level = 0.95,
  ...
)
```

## Arguments

- fit:

  A model fitted with
  [`baysesp()`](https://giovannitinervia9.github.io/baysesp/reference/baysesp.md).

- odds_ratios:

  Values of the odds ratio of verification between diseased and healthy
  subjects, \\\exp(\beta)\\; 1 is MAR. `NULL` keeps the assumption of
  the fitted model.

- unverified_priors:

  List of alternative priors on the proportion of diseased in the
  never-verified profiles. Each element is a prior, or a list with one
  prior per never-verified profile. `NULL` keeps the prior of the fitted
  model.

- rules:

  Rules, see
  [`rule_accuracy()`](https://giovannitinervia9.github.io/baysesp/reference/rule_accuracy.md).
  By default the individual tests, as in
  [`test_accuracy()`](https://giovannitinervia9.github.io/baysesp/reference/test_accuracy.md).

- level:

  Probability of the credible intervals.

- ...:

  For the structured models, arguments passed to
  [`baysesp()`](https://giovannitinervia9.github.io/baysesp/reference/baysesp.md)
  for the new fits, such as `chains` and `seed`.

## Value

A data frame with the columns `odds_ratio` (`NA` when it follows the
prior of a model with `verification = mnar()`) and `unverified_prior`,
which describe the assumptions, followed by the columns of
[`test_accuracy()`](https://giovannitinervia9.github.io/baysesp/reference/test_accuracy.md)
or
[`rule_accuracy()`](https://giovannitinervia9.github.io/baysesp/reference/rule_accuracy.md),
without `identified`.

## Details

For the structured models
([`loglinear()`](https://giovannitinervia9.github.io/baysesp/reference/loglinear.md),
[`random_effects()`](https://giovannitinervia9.github.io/baysesp/reference/random_effects.md))
the odds ratio enters the likelihood, so the model is fitted again for
each value, with the odds ratio fixed; the never-verified profiles have
no separate prior to change.

The odds ratio is fixed at each value of `odds_ratios`: with a grid of
values this is a tipping-point analysis, which shows at which odds ratio
a conclusion changes (Kosinski and Barnhart, 2003b). With
`odds_ratios = NULL` the odds ratio of the fitted model is kept (1 with
`verification = mar()`, the prior `odds_ratio_prior` with
`verification = mnar()`). The new priors in `unverified_priors` are
sampled with the random number generator of R: set the seed for
reproducible results.

**Recomputation for the saturated model**. The saturated model
estimates, for each profile \\t\\, its frequency \\\omega_t\\, the
fraction of verified subjects \\v_t\\ and the proportion of diseased
among them \\\phi^{(1)}\_t\\. With odds ratio of verification \\OR =
\exp(\beta)\\, the proportion of diseased among the unverified is
\\\phi^{(0)}\_t = \mathrm{expit}(\mathrm{logit}(\phi^{(1)}\_t) -
\beta)\\: with \\OR \> 1\\ the unverified are less often diseased than
the verified (Kosinski and Barnhart, 2003a). The proportion of diseased
in the profile is then \$\$\psi_t = v_t\\ \phi^{(1)}\_t + (1 - v_t)\\
\phi^{(0)}\_t,\$\$ or \\\psi_t = q_t\\, drawn from the prior, for a
never-verified profile. From these, \\\pi = \sum_t \omega_t \psi_t\\,
\\\theta\_{1,t} = \omega_t \psi_t / \pi\\ and \\\theta\_{0,t} = \omega_t
(1 - \psi_t) / (1 - \pi)\\, and the accuracy measures below.

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

Kosinski, A. S. and Barnhart, H. X. (2003a). Accounting for nonignorable
verification bias in assessment of diagnostic tests. *Biometrics*,
59(1), 163–171.
[doi:10.1111/1541-0420.00019](https://doi.org/10.1111/1541-0420.00019)

Kosinski, A. S. and Barnhart, H. X. (2003b). A global sensitivity
analysis of performance of a medical diagnostic test when verification
bias is present. *Statistics in Medicine*, 22(17), 2711–2721.
[doi:10.1002/sim.1517](https://doi.org/10.1002/sim.1517)

Torrance-Rynard, V. L. and Walter, S. D. (1997). Effects of dependent
errors in the assessment of diagnostic test performance. *Statistics in
Medicine*, 16(19), 2157–2175.
[doi:10.1002/(SICI)1097-0258(19971015)16:19\<2157::AID-SIM653\>3.0.CO;2-X](https://doi.org/10.1002/%28SICI%291097-0258%2819971015%2916%3A19%3C2157%3A%3AAID-SIM653%3E3.0.CO%3B2-X)

Vacek, P. M. (1985). The effect of conditional dependence on the
evaluation of diagnostic tests. *Biometrics*, 41(4), 959–968.
[doi:10.2307/2530967](https://doi.org/10.2307/2530967)

## Examples

``` r
# \donttest{
set.seed(1)
d <- simulate_study(n = 2000, prevalence = 0.2,
                    sensitivity = c(A = 0.9, B = 0.8),
                    specificity = c(A = 0.85, B = 0.95),
                    verification_prob = c(0.2, 0.9, 0.9, 1), odds_ratio = 2)
d <- baysesp_data(d, tests = c("A", "B"), verified = "verified", disease = "disease")
fit <- baysesp(d, seed = 1)
sensitivity_analysis(fit, odds_ratios = c(1, 1.5, 2, 3))
#>    odds_ratio unverified_prior   quantity      mean    median     lower
#> 1         1.0             <NA> prevalence 0.2082286 0.2079712 0.1893950
#> 2         1.0             <NA>      Se[A] 0.8762867 0.8785089 0.8260178
#> 3         1.0             <NA>      Sp[A] 0.8435602 0.8437649 0.8243615
#> 4         1.0             <NA>     PPV[A] 0.5955891 0.5955906 0.5564420
#> 5         1.0             <NA>     NPV[A] 0.9627779 0.9635954 0.9457191
#> 6         1.0             <NA>      Se[B] 0.7660146 0.7677679 0.7121700
#> 7         1.0             <NA>      Sp[B] 0.9449760 0.9451315 0.9336026
#> 8         1.0             <NA>     PPV[B] 0.7853938 0.7857336 0.7453550
#> 9         1.0             <NA>     NPV[B] 0.9387860 0.9394644 0.9215035
#> 10        1.5             <NA> prevalence 0.2048240 0.2046451 0.1864402
#> 11        1.5             <NA>      Se[A] 0.8855456 0.8872652 0.8425733
#> 12        1.5             <NA>      Sp[A] 0.8428985 0.8431162 0.8236745
#> 13        1.5             <NA>     PPV[A] 0.5921307 0.5920850 0.5529181
#> 14        1.5             <NA>     NPV[A] 0.9661582 0.9667204 0.9522381
#> 15        1.5             <NA>      Se[B] 0.7767573 0.7781269 0.7267921
#> 16        1.5             <NA>      Sp[B] 0.9447324 0.9448851 0.9334296
#> 17        1.5             <NA>     PPV[B] 0.7835142 0.7839017 0.7431078
#> 18        1.5             <NA>     NPV[B] 0.9425792 0.9430994 0.9275236
#> 19        2.0             <NA> prevalence 0.2029912 0.2028133 0.1847511
#> 20        2.0             <NA>      Se[A] 0.8904596 0.8918907 0.8518888
#> 21        2.0             <NA>      Sp[A] 0.8424890 0.8427305 0.8231370
#> 22        2.0             <NA>     PPV[A] 0.5901248 0.5901812 0.5511408
#> 23        2.0             <NA>     NPV[A] 0.9679143 0.9684844 0.9555194
#> 24        2.0             <NA>      Se[B] 0.7825642 0.7839111 0.7354772
#> 25        2.0             <NA>      Sp[B] 0.9445662 0.9447018 0.9333302
#> 26        2.0             <NA>     PPV[B] 0.7823615 0.7826878 0.7418470
#> 27        2.0             <NA>     NPV[B] 0.9445853 0.9450387 0.9303725
#> 28        3.0             <NA> prevalence 0.2010462 0.2009173 0.1832823
#> 29        3.0             <NA>      Se[A] 0.8956048 0.8969019 0.8600836
#> 30        3.0             <NA>      Sp[A] 0.8420107 0.8422476 0.8227577
#> 31        3.0             <NA>     PPV[A] 0.5878774 0.5879977 0.5488226
#> 32        3.0             <NA>     NPV[A] 0.9697255 0.9701466 0.9584984
#> 33        3.0             <NA>      Se[B] 0.7887133 0.7896355 0.7440247
#> 34        3.0             <NA>      Sp[B] 0.9443553 0.9445156 0.9331566
#> 35        3.0             <NA>     PPV[B] 0.7809997 0.7813678 0.7403787
#> 36        3.0             <NA>     NPV[B] 0.9466789 0.9471055 0.9335824
#>        upper
#> 1  0.2275395
#> 2  0.9156483
#> 3  0.8612695
#> 4  0.6342342
#> 5  0.9752787
#> 6  0.8123563
#> 7  0.9552643
#> 8  0.8227042
#> 9  0.9527518
#> 10 0.2235666
#> 11 0.9201676
#> 12 0.8605234
#> 13 0.6308545
#> 14 0.9768598
#> 15 0.8199500
#> 16 0.9551056
#> 17 0.8210365
#> 18 0.9551448
#> 19 0.2214544
#> 20 0.9225797
#> 21 0.8601487
#> 22 0.6287071
#> 23 0.9777076
#> 24 0.8241869
#> 25 0.9549732
#> 26 0.8201899
#> 27 0.9564633
#> 28 0.2193333
#> 29 0.9255052
#> 30 0.8597946
#> 31 0.6264751
#> 32 0.9787697
#> 33 0.8289792
#> 34 0.9547125
#> 35 0.8188682
#> 36 0.9579274
# }
```
