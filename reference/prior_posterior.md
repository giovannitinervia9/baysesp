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
#>  prevalence        prior 0.3827 0.3664 0.08152 0.7666
#>  prevalence    posterior 0.2192 0.2183 0.19521 0.2476
#>       Se[A]        prior 0.6464 0.7051 0.08913 0.9846
#>       Se[A]    posterior 0.8261 0.8306 0.73691 0.8905
#>       Sp[A]        prior 0.5846 0.6095 0.10014 0.9661
#>       Sp[A]    posterior 0.8388 0.8389 0.81909 0.8570
#>      PPV[A]        prior 0.5009 0.5011 0.08033 0.9181
#>      PPV[A]    posterior 0.5895 0.5898 0.54674 0.6296
#>      NPV[A]        prior 0.7377 0.7945 0.21845 0.9779
#>      NPV[A]    posterior 0.9445 0.9467 0.90768 0.9679
#>       Se[B]        prior 0.6420 0.6964 0.09622 0.9827
#>       Se[B]    posterior 0.7344 0.7373 0.65181 0.8031
#>       Sp[B]        prior 0.5845 0.6053 0.11372 0.9628
#>       Sp[B]    posterior 0.9455 0.9458 0.93359 0.9565
#>      PPV[B]        prior 0.4967 0.4896 0.07872 0.9242
#>      PPV[B]    posterior 0.7907 0.7914 0.75091 0.8302
#>      NPV[B]        prior 0.7317 0.7899 0.21307 0.9774
#>      NPV[B]    posterior 0.9265 0.9283 0.89448 0.9492
#>        q[1]        prior 0.0298 0.0269 0.00726 0.0680
#>        q[1]    posterior 0.0300 0.0273 0.00695 0.0693
plot(pp)

# }
```
