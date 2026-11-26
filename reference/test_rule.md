# Diagnostic rule class

A `test_rule` object declares each response profile positive or
negative, starting from the results of one or more tests. The accuracy
of a rule depends on the joint distribution of the tests, not only on
their sensitivities and specificities (Torrance-Rynard and Walter, 1997;
Vacek, 1985). Common rules are built with
[`single_test()`](https://giovannitinervia9.github.io/baysesp/reference/single_test.md),
[`any_positive()`](https://giovannitinervia9.github.io/baysesp/reference/any_positive.md),
[`all_positive()`](https://giovannitinervia9.github.io/baysesp/reference/all_positive.md)
and
[`min_positive()`](https://giovannitinervia9.github.io/baysesp/reference/min_positive.md);
a custom rule is built with `test_rule()`.

## Usage

``` r
test_rule(name = character(0), positive = function() NULL)
```

## Arguments

- name:

  Character string, the label of the rule in the results. Use a single
  non-missing string.

- positive:

  Function with one argument, the matrix of response profiles described
  above. It must return one logical or numeric value per row of the
  matrix, and it is applied to the matrix of the data at the moment the
  rule is used.

## Value

An object of class `test_rule`.

## Details

A response profile is the vector \\t = (t_1, \dots, t_K)\\ of the
results of the \\K\\ tests, where \\t_k = 1\\ if test \\k\\ is positive
and \\t_k = 0\\ if it is negative. With \\K\\ tests there are \\2^K\\
profiles. A rule assigns to each profile a value \\r(t)\\: \\r(t) = 1\\
if the rule declares the subject positive, \\r(t) = 0\\ otherwise.

A `test_rule` has two properties:

- `name`: the label of the rule in the results;

- `positive`: a function of the matrix of response profiles that
  returns, for each profile, whether the rule declares it positive.

The matrix of response profiles has one row per profile (\\2^K\\ rows)
and one column per test, named with the test names. Its entries are `0`
(negative) or `1` (positive). The function in `positive` must return one
value per row: `TRUE` (or `1`) if the profile is positive, `FALSE` (or
`0`) otherwise. It may refer to tests by their names, for example
`profiles[, "ct"]`.

## References

Torrance-Rynard, V. L. and Walter, S. D. (1997). Effects of dependent
errors in the assessment of diagnostic test performance. *Statistics in
Medicine*, 16(19), 2157–2175.
[doi:10.1002/(SICI)1097-0258(19971015)16:19\<2157::AID-SIM653\>3.0.CO;2-X](https://doi.org/10.1002/%28SICI%291097-0258%2819971015%2916%3A19%3C2157%3A%3AAID-SIM653%3E3.0.CO%3B2-X)

Vacek, P. M. (1985). The effect of conditional dependence on the
evaluation of diagnostic tests. *Biometrics*, 41(4), 959–968.
[doi:10.2307/2530967](https://doi.org/10.2307/2530967)

## See also

[`rule_accuracy()`](https://giovannitinervia9.github.io/baysesp/reference/rule_accuracy.md),
[`compare_tests()`](https://giovannitinervia9.github.io/baysesp/reference/compare_tests.md)

## Examples

``` r
# Positive if the tests named "ct" and "pcr" are both positive
both_ct_pcr <- test_rule(
  name = "ct and pcr positive",
  positive = function(profiles) profiles[, "ct"] == 1 & profiles[, "pcr"] == 1
)
both_ct_pcr@name
#> [1] "ct and pcr positive"
```
