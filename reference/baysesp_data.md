# Study data

Summarises a study in which `K` binary tests are performed on every
subject and the disease status is ascertained by the reference standard
only for the verified subjects (Begg and Greenes, 1983; Baker, 1995).
The data are reduced, without loss of information, to the counts of
verified diseased, verified healthy and unverified subjects for each of
the `J = 2^K` response profiles.

## Usage

``` r
baysesp_data(data, tests, verified, disease, counts = NULL)
```

## Arguments

- data:

  A data frame. By default one row per subject; with `counts`, one row
  per combination of test results, verification and disease status.

- tests:

  Names of the columns with the test results.

- verified:

  Name of the column that says whether the subject was verified.

- disease:

  Name of the column with the disease status. It may be missing (`NA`)
  only for unverified subjects, and it is ignored for them.

- counts:

  Optional name of a column with the number of subjects in each row, for
  data that are already aggregated.

## Value

An object of class `baysesp_data`.

## Details

Test results, verification and disease status must be coded as `0`/`1`
or as logical values, with `1` (`TRUE`) meaning positive, verified and
diseased.

The object has the properties `@tests`, `@profiles` (a `J x K` matrix of
the response profiles, in binary order: the first row has all tests
negative and the first test varies fastest), `@n_diseased`, `@n_healthy`
and `@n_unverified` (counts for each profile) and the read-only
`@never_verified` (profiles without any verified subject).
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the aggregated table.

## References

Baker, S. G. (1995). Evaluating multiple diagnostic tests with partial
verification. *Biometrics*, 51(1), 330–337.
[doi:10.2307/2533339](https://doi.org/10.2307/2533339)

Begg, C. B. and Greenes, R. A. (1983). Assessment of diagnostic tests
when disease verification is subject to selection bias. *Biometrics*,
39(1), 207–215. [doi:10.2307/2530820](https://doi.org/10.2307/2530820)

## Examples

``` r
d <- simulate_study(n = 1000, prevalence = 0.1,
                    sensitivity = c(A = 0.9, B = 0.8),
                    specificity = c(A = 0.95, B = 0.9),
                    verification_prob = c(0.1, 0.8, 0.8, 1))
baysesp_data(d, tests = c("A", "B"), verified = "verified", disease = "disease")
#> <baysesp_data> 2 tests, 4 profiles, 1000 subjects, 291 verified
#> never-verified profiles: none
#> 
#>  A B diseased healthy unverified
#>  0 0        0      88        683
#>  1 0       14      30         12
#>  0 1        7      74         14
#>  1 1       72       6          0

# aggregated data
agg <- data.frame(T1 = c(0, 1, 0, 1, 0, 1),
                  V  = c(1, 1, 1, 1, 0, 0),
                  D  = c(0, 0, 1, 1, NA, NA),
                  n  = c(40, 5, 2, 30, 300, 20))
baysesp_data(agg, tests = "T1", verified = "V", disease = "D", counts = "n")
#> <baysesp_data> 1 test, 2 profiles, 397 subjects, 77 verified
#> never-verified profiles: none
#> 
#>  T1 diseased healthy unverified
#>   0        2      40        300
#>   1       30       5         20
```
