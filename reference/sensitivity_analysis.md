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
#> 1         1.0             <NA> prevalence 0.2092226 0.2089931 0.1902462
#> 2         1.0             <NA>      Se[A] 0.8712684 0.8733881 0.8215693
#> 3         1.0             <NA>      Sp[A] 0.8430957 0.8431539 0.8248153
#> 4         1.0             <NA>     PPV[A] 0.5949234 0.5948473 0.5563297
#> 5         1.0             <NA>     NPV[A] 0.9610812 0.9618838 0.9437379
#> 6         1.0             <NA>      Se[B] 0.7617043 0.7626477 0.7083632
#> 7         1.0             <NA>      Sp[B] 0.9447449 0.9449223 0.9332584
#> 8         1.0             <NA>     PPV[B] 0.7847668 0.7851835 0.7443982
#> 9         1.0             <NA>     NPV[B] 0.9373738 0.9378300 0.9205118
#> 10        1.5             <NA> prevalence 0.2055242 0.2054128 0.1875115
#> 11        1.5             <NA>      Se[A] 0.8816344 0.8831522 0.8389579
#> 12        1.5             <NA>      Sp[A] 0.8424899 0.8425014 0.8242849
#> 13        1.5             <NA>     PPV[A] 0.5914526 0.5913568 0.5529788
#> 14        1.5             <NA>     NPV[A] 0.9648806 0.9653994 0.9505215
#> 15        1.5             <NA>      Se[B] 0.7734285 0.7739942 0.7246175
#> 16        1.5             <NA>      Sp[B] 0.9445232 0.9446663 0.9330735
#> 17        1.5             <NA>     PPV[B] 0.7828895 0.7833079 0.7423628
#> 18        1.5             <NA>     NPV[B] 0.9415361 0.9417545 0.9268275
#> 19        2.0             <NA> prevalence 0.2035421 0.2034370 0.1859001
#> 20        2.0             <NA>      Se[A] 0.8871343 0.8882654 0.8478650
#> 21        2.0             <NA>      Sp[A] 0.8421083 0.8421565 0.8240588
#> 22        2.0             <NA>     PPV[A] 0.5894395 0.5892667 0.5509857
#> 23        2.0             <NA>     NPV[A] 0.9668492 0.9672430 0.9545814
#> 24        2.0             <NA>      Se[B] 0.7797582 0.7802408 0.7331961
#> 25        2.0             <NA>      Sp[B] 0.9443679 0.9444936 0.9329572
#> 26        2.0             <NA>     PPV[B] 0.7817375 0.7821610 0.7413572
#> 27        2.0             <NA>     NPV[B] 0.9437296 0.9438987 0.9302430
#> 28        3.0             <NA> prevalence 0.2014460 0.2013237 0.1840445
#> 29        3.0             <NA>      Se[A] 0.8928895 0.8937604 0.8573904
#> 30        3.0             <NA>      Sp[A] 0.8416580 0.8417122 0.8236370
#> 31        3.0             <NA>     PPV[A] 0.5871837 0.5870358 0.5489598
#> 32        3.0             <NA>     NPV[A] 0.9688750 0.9692224 0.9579226
#> 33        3.0             <NA>      Se[B] 0.7864543 0.7867841 0.7424804
#> 34        3.0             <NA>      Sp[B] 0.9441679 0.9442843 0.9327506
#> 35        3.0             <NA>     PPV[B] 0.7803753 0.7807117 0.7397278
#> 36        3.0             <NA>     NPV[B] 0.9460126 0.9461758 0.9333552
#>        upper
#> 1  0.2288867
#> 2  0.9108137
#> 3  0.8598069
#> 4  0.6333556
#> 5  0.9737539
#> 6  0.8085791
#> 7  0.9554404
#> 8  0.8233076
#> 9  0.9515504
#> 10 0.2239159
#> 11 0.9163834
#> 12 0.8592548
#> 13 0.6300658
#> 14 0.9757131
#> 15 0.8177608
#> 16 0.9552004
#> 17 0.8213861
#> 18 0.9541086
#> 19 0.2216201
#> 20 0.9194891
#> 21 0.8588371
#> 22 0.6280094
#> 23 0.9768072
#> 24 0.8227492
#> 25 0.9550404
#> 26 0.8199734
#> 27 0.9555378
#> 28 0.2193324
#> 29 0.9229309
#> 30 0.8583656
#> 31 0.6255749
#> 32 0.9779935
#> 33 0.8277532
#> 34 0.9548433
#> 35 0.8185556
#> 36 0.9570919
# }
```
