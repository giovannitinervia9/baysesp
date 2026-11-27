# Compare two tests or rules

Posterior summaries of the differences in sensitivity, specificity and
Youden index (`first` minus `second`) and of the ratios of the true and
false positive rates (`first` over `second`) (Alonzo, 2005; Pepe and
Alonzo, 2001; Schatzkin et al., 1987).

## Usage

``` r
compare_tests(fit, first, second, level = 0.95)
```

## Arguments

- fit:

  A fitted model, from
  [`baysesp()`](https://giovannitinervia9.github.io/baysesp/reference/baysesp.md).

- first, second:

  Two tests (by name) or two rules, see
  [test_rule](https://giovannitinervia9.github.io/baysesp/reference/test_rule.md).

- level:

  Probability of the credible intervals.

## Value

A data frame as in
[`test_accuracy()`](https://giovannitinervia9.github.io/baysesp/reference/test_accuracy.md),
with the additional column `probability`: the posterior probability that
a difference is positive or that a ratio is greater than 1.

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

For two rules \\R\\ (`first`) and \\R'\\ (`second`):

- **differences**: \\\Delta Se = Se_R - Se\_{R'}\\, \\\Delta Sp = Sp_R -
  Sp\_{R'}\\ and, for the Youden index \\J_R = Se_R + Sp_R - 1\\, which
  summarizes the accuracy in one number (0 for a useless test, 1 for a
  perfect one), \\\Delta J = \Delta Se + \Delta Sp\\. A positive
  difference favors the first rule;

- **ratios** of the true positive rates, \\rTPR = Se_R / Se\_{R'}\\, and
  of the false positive rates, \\rFPR = (1 - Sp_R) / (1 - Sp\_{R'})\\.
  For example, \\rTPR = 1.2\\ means that the first rule detects 20% more
  diseased subjects; \\rFPR \< 1\\ means that it gives fewer false
  positives.

The column `probability` is the posterior probability that a difference
is positive, or that a ratio is greater than 1. \\P(rTPR \> 1)\\ equals
\\P(\Delta Se \> 0)\\, and \\P(rFPR \> 1)\\ equals \\P(\Delta Sp \<
0)\\.

When some profile is never verified, the differences depend on the prior
of its proportion of diseased (Gustafson, 2015), but the ratios do not,
as long as neither rule declares a never-verified profile positive
(Alonzo, 2005; Pepe and Alonzo, 2001; Schatzkin et al., 1987). In that
case the signs of the differences in sensitivity and in specificity are
also identified, through the probabilities of the ratios.

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
compare_tests(fit, "A", "B")
#>                quantity          mean        median       lower       upper
#> 1         Se[A] - Se[B]  0.0973789883  0.0964292438  0.02417792  0.17155312
#> 2         Sp[A] - Sp[B] -0.0975348395 -0.0975532413 -0.12736462 -0.06761717
#> 3 Youden[A] - Youden[B] -0.0001558511 -0.0003221472 -0.08388671  0.08097892
#> 4       TPR[A] / TPR[B]  1.1237196494  1.1211319813  1.02874656  1.23025364
#> 5       FPR[A] / FPR[B]  2.9172331397  2.8595559824  2.04731457  4.08945916
#>   identified probability
#> 1       TRUE     0.99625
#> 2       TRUE     0.00000
#> 3       TRUE     0.49775
#> 4       TRUE     0.99625
#> 5       TRUE     1.00000
compare_tests(fit, any_positive(), "A")
#>                           quantity        mean      median        lower
#> 1         Se[any positive] - Se[A]  0.09420604  0.09248344  0.056695901
#> 2         Sp[any positive] - Sp[A] -0.04598212 -0.04553350 -0.061961373
#> 3 Youden[any positive] - Youden[A]  0.04822392  0.04719800  0.006315887
#> 4       TPR[any positive] / TPR[A]  1.10589169  1.10319126  1.061246362
#> 5       FPR[any positive] / FPR[A]  1.30868275  1.30382676  1.206844399
#>         upper identified probability
#> 1  0.13898661       TRUE       1.000
#> 2 -0.03252576       TRUE       0.000
#> 3  0.09598884       TRUE       0.991
#> 4  1.16270481       TRUE       1.000
#> 5  1.43781302       TRUE       1.000
# }
```
