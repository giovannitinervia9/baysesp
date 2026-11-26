# Missing at random verification

`mar()` is the assumption that verification depends only on the test
results. Once the response profile is known, the disease status does not
change the probability of being verified. The function has no arguments.

## Usage

``` r
mar()
```

## Value

An object of class `mar`.

## Details

Let \\V = 1\\ indicate that a subject is verified, \\T\\ be the response
profile and \\D\\ the disease status, and let \\\lambda\_{t,d} = P(V = 1
\mid T = t, D = d)\\. The assumption is \$\$\lambda\_{t,1} =
\lambda\_{t,0} = \lambda_t \quad \textrm{for every profile } t,\$\$ so
the probability of verification depends on the profile only.

Equivalently, among subjects with the same profile \\t\\, verified and
unverified subjects have the same proportion of diseased: \$\$P(D = 1
\mid T = t, V = 1) = P(D = 1 \mid T = t, V = 0) = P(D = 1 \mid T =
t).\$\$ Verification therefore carries no information on the disease
beyond the test results.

Under this assumption the verification process is ignorable and can be
left out of the likelihood. Let \\n\_{t,d}\\ be the number of verified
subjects with profile \\t\\ and disease status \\d\\, \\u_t\\ the number
of unverified subjects with profile \\t\\, \\\pi = P(D = 1)\\ the
prevalence, \\\pi_1 = \pi\\, \\\pi_0 = 1 - \pi\\, and \\\theta\_{d,t} =
P(T = t \mid D = d)\\. The log-likelihood of the observed counts is
\$\$\ell = \sum_t \Big\[\sum_d n\_{t,d}
\log(\pi_d\\\theta\_{d,t}\\\lambda\_{t,d}) + u_t \log \sum_d
\pi_d\\\theta\_{d,t}\\(1 - \lambda\_{t,d})\Big\].\$\$ Setting
\\\lambda\_{t,d} = \lambda_t\\ splits it into two parts, \$\$\ell =
\ell_D + \ell_V,\$\$ \$\$\ell_D = \sum_t \Big\[\sum_d n\_{t,d}
\log(\pi_d\\\theta\_{d,t}) + u_t \log \sum_d \pi_d\\\theta\_{d,t}\Big\],
\qquad \ell_V = \sum_t \Big\[(n\_{t,1} + n\_{t,0}) \log \lambda_t + u_t
\log(1 - \lambda_t)\Big\].\$\$ The verification probabilities
\\\lambda_t\\ appear only in \\\ell_V\\. Inference on the prevalence and
on the test distributions is based on \\\ell_D\\, and the unverified
subjects contribute to it through the totals \\u_t\\.

MAR is plausible when the choice of who undergoes the reference standard
is made from the test results alone, for example when a protocol
prescribes the reference test for every subject with a given test
pattern, and no other information related to the disease affects the
choice. It is doubtful when the decision to verify uses information on
the disease that the tests do not capture, such as symptoms, medical
history or clinical suspicion, or when, within the same profile,
subjects who are more likely to be diseased are verified more often. In
those cases the proportion of diseased among verified subjects differs
from that among unverified subjects with the same profile, and the
estimates obtained under MAR are biased. The data do not contain
information on whether MAR holds, so its plausibility must be judged on
substantive grounds. The alternative in which verification also depends
on the disease status is missing not at random.

Profiles that are never verified are compatible with MAR, but for them
the data give no information on the proportion of diseased: there are no
verified subjects from which to estimate it. What is concluded about it
comes entirely from an assumption: a prior on that proportion, or a
structure imposed on the dependence between the tests.

## See also

[`mnar()`](https://giovannitinervia9.github.io/baysesp/reference/mnar.md),
[verification_mechanism](https://giovannitinervia9.github.io/baysesp/reference/verification_mechanism.md)

## Examples

``` r
mar()
#> <mar> verification depends only on the test results
```
