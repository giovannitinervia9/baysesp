# Single test rule

A rule that declares a subject positive when one given test is positive,
whatever the results of the other tests.

## Usage

``` r
single_test(test)
```

## Arguments

- test:

  Character string, the name of one test. It must be a single
  non-missing string that matches a test name in the data. There is no
  default.

## Value

An object of class `test_rule`.

## Details

A response profile is the vector \\t = (t_1, \dots, t_K)\\ of the
results of the \\K\\ tests, where \\t_k = 1\\ if test \\k\\ is positive
and \\t_k = 0\\ if it is negative. With \\K\\ tests there are \\2^K\\
profiles.

`single_test(test)` gives the rule \$\$r(t) = t_k,\$\$ where \\k\\ is
the test named `test`. The rule declares positive the profiles in which
that test is positive, and ignores the results of the other tests. Its
label in the results is the name of the test.

Functions that accept rules also accept a single test name in place of a
rule: the name is then used as a rule of this kind.

The test name is checked when the rule is used, not when it is created:
if the data do not contain a test with that name, an error is raised.

## See also

[`rule_accuracy()`](https://giovannitinervia9.github.io/baysesp/reference/rule_accuracy.md),
[`compare_tests()`](https://giovannitinervia9.github.io/baysesp/reference/compare_tests.md)

## Examples

``` r
single_test("ct")
#> <test_rule> ct
```
