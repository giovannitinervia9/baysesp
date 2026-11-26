# Bayesian estimation of the accuracy of diagnostic tests

Fits a Bayesian model for the sensitivity and specificity of one or more
binary tests when not all subjects have been verified with the reference
standard. The objects of the model classes passed to `model` contain the
specification of the priors and of the verification mechanism.

## Usage

``` r
baysesp(
  data,
  model = saturated(),
  backend = c("rstan", "cmdstanr"),
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = NULL,
  ...
)
```

## Arguments

- data:

  A `baysesp_data` object, created by
  [`baysesp_data()`](https://giovannitinervia9.github.io/baysesp/reference/baysesp_data.md),
  with the results of the tests, the verification indicator and the
  disease status of each subject. Any other object gives an error.

- model:

  An accuracy model: an object created by one of the model constructors
  [`saturated()`](https://giovannitinervia9.github.io/baysesp/reference/saturated.md),
  [`loglinear()`](https://giovannitinervia9.github.io/baysesp/reference/loglinear.md)
  or
  [`random_effects()`](https://giovannitinervia9.github.io/baysesp/reference/random_effects.md).
  It holds the priors on the model parameters and the verification
  mechanism, which is `verification = mar()` or `verification = mnar()`
  given to the constructor. Default
  [`saturated()`](https://giovannitinervia9.github.io/baysesp/reference/saturated.md),
  with MAR verification. Any other object gives an error.

- backend:

  Stan interface, a single value: `"rstan"` (default) or `"cmdstanr"`.
  With `"rstan"`, the Stan programs are the ones compiled when the
  package is installed, so nothing is compiled at run time. With
  `"cmdstanr"`, each program is compiled the first time it is used and
  the executable is kept in the user cache directory, so the first run
  is slower and later runs reuse the compiled program. If `cmdstanr` or
  CmdStan is not available, `"rstan"` is used and a message is printed.
  The `backend` element of the result reports the backend that was
  actually used. Any other value gives an error.

- chains:

  Positive integer, the number of Markov chains, which run
  independently. Default `4`.

- iter_warmup:

  Positive integer, the number of warmup iterations of each chain. The
  warmup iterations adapt the sampler and are then discarded. Default
  `1000`.

- iter_sampling:

  Positive integer, the number of sampling iterations of each chain,
  kept as draws. Default `1000`, which gives `chains * 1000` draws in
  total.

- seed:

  `NULL` (default) or a single whole number between `0` and
  `2147483647`. If `NULL`, a seed is drawn from the random number
  generator of R, so that calling
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before `baysesp()`
  makes the run reproducible. A number is used directly as the seed of
  the sampler.

- ...:

  Further arguments passed to
  [`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html)
  when `backend` is `"rstan"` (for example `cores`, to run the chains in
  parallel, or `control`, for example `list(adapt_delta = 0.95)`), or to
  the `$sample()` method of `cmdstanr` when `backend` is `"cmdstanr"`
  (for example `parallel_chains`, or `adapt_delta`). Do not repeat
  `chains`, `iter_warmup`, `iter_sampling` or `seed`, which `baysesp()`
  already passes to the backend.

## Value

An object of class `baysesp_fit` with the following properties:

- `model`: the model object given in `model`.

- `data`: the `baysesp_data` object given in `data`.

- `draws`: a
  [`posterior::draws_array`](https://mc-stan.org/posterior/reference/draws_array.html)
  with the sampling iterations only, that is `iter_sampling` iterations,
  `chains` chains and the variables listed below.

- `backend`: the backend used, `"rstan"` or `"cmdstanr"`.

- `stan_fit`: the object returned by the backend: a `stanfit` object for
  `"rstan"`, a `CmdStanMCMC` object for `"cmdstanr"`.

The profiles are the \\J = 2^K\\ response profiles of the \\K\\ tests,
indexed in binary order: profile \\j\\ is the pattern of results with
index \\j = 1 + \sum_k t_k 2^{k-1}\\, where \\t_k = 1\\ for a positive
result on test \\k\\ and \\t_k = 0\\ otherwise. Profile 1 has all tests
negative, and the first test changes fastest. The variables in `draws`
are:

- `prevalence`: the prevalence of the disease.

- `theta_healthy[j]` and `theta_diseased[j]`: the probability of profile
  \\j\\ among the healthy and among the diseased. Within each draw, each
  of the two sets sums to 1 over \\j\\.

- `omega[j]`, `v[j]` and `phi1[j]`, only for the saturated model: the
  probability of profile \\j\\, the probability of verification given
  profile \\j\\, and the probability of disease among the verified
  subjects with profile \\j\\. `v[j]` and `phi1[j]` are 0 for profiles
  with no verified subjects.

- `q[b]`: the prevalence of disease among the subjects with the \\b\\-th
  profile that has no verified subjects. Present only if some profile
  has no verified subjects.

- `beta`: the log odds ratio of verification between diseased and
  healthy subjects with the same profile, present only with MNAR
  verification. Without `by_profile = TRUE` it is a single value. With
  `by_profile = TRUE` there is one value, named `beta[i]`, for each
  profile with verified subjects, in binary order.

## Details

The posterior distribution of the parameters of `model`, given `data`,
is approximated with the No-U-Turn sampler of Stan, a form of
Hamiltonian Monte Carlo that adapts its step size and mass matrix during
the warmup iterations. The likelihood and the priors are coded in a Stan
program specific to the model class and to the verification mechanism.
Each chain runs `iter_warmup` warmup iterations, which are discarded,
followed by `iter_sampling` sampling iterations, which are kept as
draws. The draws therefore form `chains * iter_sampling` joint samples
of the parameters of the model, and every quantity reported by the
package is computed from them.

## See also

[`test_accuracy()`](https://giovannitinervia9.github.io/baysesp/reference/test_accuracy.md)

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
test_accuracy(fit)
#>     quantity      mean    median     lower     upper identified
#> 1 prevalence 0.2052042 0.2048145 0.1795179 0.2316453       TRUE
#> 2      Se[A] 0.8851116 0.8871800 0.8226758 0.9317420       TRUE
#> 3      Sp[A] 0.8488080 0.8492479 0.8232210 0.8731027       TRUE
#> 4     PPV[A] 0.6017400 0.6018114 0.5445422 0.6591683       TRUE
#> 5     NPV[A] 0.9661723 0.9671865 0.9458064 0.9805380       TRUE
#> 6      Se[B] 0.7878260 0.7892088 0.7192221 0.8507625       TRUE
#> 7      Sp[B] 0.9465961 0.9470981 0.9301900 0.9614485       TRUE
#> 8     PPV[B] 0.7920101 0.7926048 0.7341592 0.8445836       TRUE
#> 9     NPV[B] 0.9452112 0.9460813 0.9240182 0.9626702       TRUE
# }
```
