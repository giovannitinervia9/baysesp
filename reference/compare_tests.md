# Compare two tests or rules

Posterior summaries of the differences in sensitivity, specificity and
Youden index (`first` minus `second`) and of the ratios of the true and
false positive rates (`first` over `second`).

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
their individual sensitivities and specificities.

Each quantity is computed from every posterior draw of \\\pi\\,
\\\theta_0\\ and \\\theta_1\\. The results report its posterior mean and
median and an equal-tailed credible interval with probability `level`.

For two rules \\R\\ (`first`) and \\R'\\ (`second`):

- **differences**: \\\Delta Se = Se_R - Se\_{R'}\\, \\\Delta Sp = Sp_R -
  Sp\_{R'}\\ and, for the Youden index \\J_R = Se_R + Sp_R - 1\\, which
  summarises the accuracy in one number (0 for a useless test, 1 for a
  perfect one), \\\Delta J = \Delta Se + \Delta Sp\\. A positive
  difference favours the first rule;

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
of its proportion of diseased, but the ratios do not, as long as neither
rule declares a never-verified profile positive. In that case the signs
of the differences in sensitivity and in specificity are also
identified, through the probabilities of the ratios.

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
#> 1         Se[A] - Se[B]  0.0972856186  0.0977011712  0.02168599  0.17650180
#> 2         Sp[A] - Sp[B] -0.0977881666 -0.0975503042 -0.12733509 -0.06859078
#> 3 Youden[A] - Youden[B] -0.0005025481 -0.0008030661 -0.08393341  0.08597781
#> 4       TPR[A] / TPR[B]  1.1252383139  1.1232906839  1.02670754  1.23966427
#> 5       FPR[A] / FPR[B]  2.8959639962  2.8492295512  2.05311774  4.02143271
#>   identified probability
#> 1       TRUE      0.9920
#> 2       TRUE      0.0000
#> 3       TRUE      0.4945
#> 4       TRUE      0.9920
#> 5       TRUE      1.0000
compare_tests(fit, any_positive(), "A")
#>                           quantity        mean      median        lower
#> 1         Se[any positive] - Se[A]  0.09423374  0.09274048  0.055726849
#> 2         Sp[any positive] - Sp[A] -0.04595762 -0.04552620 -0.061781377
#> 3 Youden[any positive] - Youden[A]  0.04827612  0.04666800  0.006056152
#> 4       TPR[any positive] / TPR[A]  1.10710037  1.10470621  1.060645167
#> 5       FPR[any positive] / FPR[A]  1.30643854  1.30173813  1.205143968
#>         upper identified probability
#> 1  0.14027900       TRUE     1.00000
#> 2 -0.03222708       TRUE     0.00000
#> 3  0.09702180       TRUE     0.98775
#> 4  1.16696488       TRUE     1.00000
#> 5  1.43132449       TRUE     1.00000
# }
```
