# Prior distributions

`prior_distribution` is the abstract parent class of all prior
distributions. Its concrete subclasses are the Beta, normal, Student t
and Dirichlet priors, each with its own constructor.

## Usage

``` r
prior_distribution(family = character(0), support = character(0))
```

## Arguments

- family:

  Name of the distribution family (read-only property).

- support:

  Scale the prior acts on (read-only property): `"probability"` for
  priors on a probability, `"real"` for priors on an unconstrained
  parameter, `"simplex"` for priors on a vector of probabilities that
  sum to one.

## See also

[`prior_beta()`](https://giovannitinervia9.github.io/baysesp/reference/prior_beta.md),
[`elicit_beta()`](https://giovannitinervia9.github.io/baysesp/reference/elicit_beta.md),
[`prior_normal()`](https://giovannitinervia9.github.io/baysesp/reference/prior_normal.md),
[`prior_student()`](https://giovannitinervia9.github.io/baysesp/reference/prior_student.md),
[`prior_dirichlet()`](https://giovannitinervia9.github.io/baysesp/reference/prior_dirichlet.md)
