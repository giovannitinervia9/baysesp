# All positive rule (AND rule)

A rule that declares a subject positive only when all the tests are
positive.

## Usage

``` r
all_positive()
```

## Value

An object of class `test_rule`.

## Details

A response profile is the vector \\t = (t_1, \dots, t_K)\\ of the
results of the \\K\\ tests, where \\t_k = 1\\ if test \\k\\ is positive
and \\t_k = 0\\ if it is negative. With \\K\\ tests there are \\2^K\\
profiles.

`all_positive()` gives the rule \$\$r(t) = 1 \textrm{ if }
\sum\_{k=1}^{K} t_k = K, \quad r(t) = 0 \textrm{ otherwise.}\$\$ The
only profile declared positive is the one in which all the tests are
positive, \\t = (1, \dots, 1)\\. The label of the rule in the results is
`"all positive"`.

This function has no arguments.

## See also

[`rule_accuracy()`](https://giovannitinervia9.github.io/baysesp/reference/rule_accuracy.md),
[`compare_tests()`](https://giovannitinervia9.github.io/baysesp/reference/compare_tests.md)

## Examples

``` r
all_positive()
#> <test_rule> all positive
```
