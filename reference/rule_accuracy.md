# Accuracy of diagnostic rules

Posterior summaries of the sensitivity (`Se`), specificity (`Sp`),
positive predictive value (`PPV`) and negative predictive value (`NPV`)
of rules that combine the tests. For a single test the result is the
same as in
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
their individual sensitivities and specificities.

Each quantity is computed from every posterior draw of \\\pi\\,
\\\theta_0\\ and \\\theta_1\\. The results report its posterior mean and
median and an equal-tailed credible interval with probability `level`.

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
#> 1  Se[any positive] 0.9793453 0.9851305 0.9290624 0.9994905       TRUE
#> 2  Sp[any positive] 0.8028503 0.8032343 0.7741721 0.8298526       TRUE
#> 3 PPV[any positive] 0.5618369 0.5622735 0.5090977 0.6144496       TRUE
#> 4 NPV[any positive] 0.9933103 0.9952834 0.9758382 0.9998439       TRUE
#> 5  Se[all positive] 0.6935923 0.6944468 0.6214782 0.7582623       TRUE
#> 6  Sp[all positive] 0.9925537 0.9929214 0.9858199 0.9972360       TRUE
#> 7 PPV[all positive] 0.9600218 0.9620041 0.9251793 0.9849610       TRUE
#> 8 NPV[all positive] 0.9261159 0.9267453 0.9048655 0.9443424       TRUE
# }
```
