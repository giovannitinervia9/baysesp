# Compare prior and posterior

Fits the same model again sampling from the prior alone, and compares
the prior and posterior distributions of the accuracy of the tests (or
of rules) and, for the saturated model, of the parameters the data
cannot estimate: the proportion of diseased in never-verified profiles
(`q`) and the log odds ratio of verification (`beta`) when
`verification = mnar()` (Kosinski and Barnhart, 2003). This shows how
much each result depends on the data and how much on the assumptions: a
posterior equal to its prior is not learned from the data (Buzoianu and
Kadane, 2008; Gustafson, 2015). For the structured models it is also a
prior simulation, which shows the prior that their priors induce on the
accuracy.

## Usage

``` r
prior_posterior(fit, rules = NULL, level = 0.95, ...)
```

## Arguments

- fit:

  A model fitted with
  [`baysesp()`](https://giovannitinervia9.github.io/baysesp/reference/baysesp.md).

- rules:

  Rules, see
  [`rule_accuracy()`](https://giovannitinervia9.github.io/baysesp/reference/rule_accuracy.md).
  By default the individual tests, as in
  [`test_accuracy()`](https://giovannitinervia9.github.io/baysesp/reference/test_accuracy.md).

- level:

  Probability of the credible intervals.

- ...:

  Arguments passed to
  [`baysesp()`](https://giovannitinervia9.github.io/baysesp/reference/baysesp.md)
  for the prior fit, such as `chains` and `seed`.

## Value

A data frame of class `prior_posterior` with columns `quantity`,
`distribution` (`"prior"` or `"posterior"`), `mean`, `median`, `lower`
and `upper`. It has a
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) method.

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

Buzoianu, M. and Kadane, J. B. (2008). Adjusting for verification bias
in diagnostic test evaluation: A Bayesian approach. *Statistics in
Medicine*, 27(13), 2453–2473.
[doi:10.1002/sim.3099](https://doi.org/10.1002/sim.3099)

Gustafson, P. (2015). *Bayesian Inference for Partially Identified
Models: Exploring the Limits of Limited Data*. Chapman & Hall/CRC, Boca
Raton, FL.

Kosinski, A. S. and Barnhart, H. X. (2003). Accounting for nonignorable
verification bias in assessment of diagnostic tests. *Biometrics*,
59(1), 163–171.
[doi:10.1111/1541-0420.00019](https://doi.org/10.1111/1541-0420.00019)

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
                    verification_prob = c(0, 0.9, 0.9, 1))
d <- baysesp_data(d, tests = c("A", "B"), verified = "verified", disease = "disease")
fit <- baysesp(d, saturated(unverified_prevalence_prior = elicit_beta(0.03, 0.06)),
               seed = 1)
pp <- prior_posterior(fit, seed = 1)
pp
#>    quantity distribution   mean median   lower  upper
#>  prevalence        prior 0.3769 0.3525 0.05291 0.8161
#>  prevalence    posterior 0.2191 0.2184 0.19574 0.2469
#>       Se[A]        prior 0.6369 0.7104 0.04289 0.9901
#>       Se[A]    posterior 0.8256 0.8293 0.73899 0.8910
#>       Sp[A]        prior 0.5977 0.6245 0.08895 0.9856
#>       Sp[A]    posterior 0.8388 0.8391 0.81930 0.8572
#>      PPV[A]        prior 0.4968 0.4865 0.02429 0.9764
#>      PPV[A]    posterior 0.5893 0.5894 0.54859 0.6287
#>      NPV[A]        prior 0.7397 0.8227 0.16199 0.9890
#>      NPV[A]    posterior 0.9444 0.9463 0.90977 0.9678
#>       Se[B]        prior 0.6319 0.7031 0.04945 0.9917
#>       Se[B]    posterior 0.7340 0.7362 0.65379 0.8005
#>       Sp[B]        prior 0.5989 0.6344 0.08157 0.9839
#>       Sp[B]    posterior 0.9457 0.9461 0.93364 0.9570
#>      PPV[B]        prior 0.4937 0.4889 0.02383 0.9785
#>      PPV[B]    posterior 0.7911 0.7917 0.74845 0.8317
#>      NPV[B]        prior 0.7373 0.8200 0.15561 0.9901
#>      NPV[B]    posterior 0.9265 0.9277 0.89612 0.9484
#>        q[1]        prior 0.0303 0.0276 0.00742 0.0693
#>        q[1]    posterior 0.0302 0.0278 0.00765 0.0671
plot(pp)

# }
```
