# Dirichlet prior

A Dirichlet distribution for a vector of probabilities that sum to one,
such as the frequencies of the test result profiles. It generalizes the
Beta distribution to more than two components: each component is Beta
distributed on its own.

## Usage

``` r
prior_dirichlet(concentration = 1)
```

## Arguments

- concentration:

  Concentration parameters, positive finite numbers. Either a single
  number, repeated for every component, or a numeric vector with one
  value per component (one per response profile, so its length is 1 or
  \\J = 2^K\\, where \\K\\ is the number of tests). The length is
  checked when the prior is given to a model. Default 1, the uniform
  distribution on the simplex.

## Value

An object of class `prior_dirichlet`.

## Details

Let \\a = (a_1, \dots, a_J)\\ be the concentration parameters and \\p =
(p_1, \dots, p_J)\\, with \\\sum_j p_j = 1\\, the vector of
probabilities. The distribution has the following properties:

- the prior mean of component \\j\\ is \\a_j / \sum_k a_k\\;

- the total \\\sum_k a_k\\ acts as a prior sample size: the larger it
  is, the more concentrated the distribution is around its mean, while a
  small total leaves the probabilities widely spread;

- each component has marginal distribution \\\mathrm{Beta}(a_j, \sum_k
  a_k - a_j)\\.

A single value \\c\\ repeated for all components gives a symmetric prior
with mean \\1/J\\ for every component. With all concentrations equal to
1 the distribution is uniform on the simplex of probabilities.
Concentrations below 1 move the mass towards the corners of the simplex,
where some components are close to zero; larger concentrations move it
towards the prior mean.

## Examples

``` r
prior_dirichlet(1)
#> <prior_dirichlet> Dirichlet(concentration = 1)
#>   acts on a vector of probabilities
prior_dirichlet(c(2, 1, 1, 1))
#> <prior_dirichlet> Dirichlet(concentration = 2, 1, 1, 1)
#>   acts on a vector of probabilities
```
