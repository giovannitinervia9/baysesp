# Plot the posterior predictive checks of a model comparison

For each model, the observed count of each cell of the table (verified
diseased, verified healthy and unverified subjects of each response
profile of the tests) is shown with the interval of its replicated
counts under the posterior predictive distribution. Observed counts
outside their intervals are cells that the model does not reproduce. The
counts are on a square root scale, so that small and large cells are
both visible. The label of each model reports the posterior predictive
p-value of the chi-square discrepancy of the whole table.

## Usage

``` r
# S3 method for class 'model_comparison'
plot(x, ...)
```

## Arguments

- x:

  Result of
  [`compare_models()`](https://giovannitinervia9.github.io/baysesp/reference/compare_models.md).

- ...:

  Not used.

## Value

A `ggplot` object.
