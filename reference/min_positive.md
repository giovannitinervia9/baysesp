# Minimum positive rule

A rule that declares a subject positive when at least `k` of the tests
are positive.

## Usage

``` r
min_positive(k)
```

## Arguments

- k:

  Single whole number, at least 1, given as a number (for example `2`).
  It must not exceed the number of tests in the data: this is checked
  when the rule is used, and a larger value raises an error. There is no
  default.

## Value

An object of class `test_rule`.

## Details

A response profile is the vector \\t = (t_1, \dots, t_K)\\ of the
results of the \\K\\ tests, where \\t_k = 1\\ if test \\k\\ is positive
and \\t_k = 0\\ if it is negative. With \\K\\ tests there are \\2^K\\
profiles.

`min_positive(k)` gives the rule \$\$r(t) = 1 \textrm{ if }
\sum\_{k=1}^{K} t_k \ge k, \quad r(t) = 0 \textrm{ otherwise.}\$\$ With
`k = 1` the rule declares positive every profile with at least one
positive test; with `k` equal to the number of tests it declares
positive only the profile with all the tests positive. The label of the
rule in the results is `"at least k positive"`, with the value of `k`.

## See also

[`rule_accuracy()`](https://giovannitinervia9.github.io/baysesp/reference/rule_accuracy.md),
[`compare_tests()`](https://giovannitinervia9.github.io/baysesp/reference/compare_tests.md)

## Examples

``` r
min_positive(2)
#> <test_rule> at least 2 positive
```
