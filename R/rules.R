#' @include generics.R
NULL

#' Diagnostic rule class
#'
#' A `test_rule` object declares each response profile positive or negative,
#' starting from the results of one or more tests. The accuracy of a rule
#' depends on the joint distribution of the tests, not only on their
#' sensitivities and specificities. Common rules are built with `single_test()`,
#' `any_positive()`, `all_positive()` and `min_positive()`; a custom rule is
#' built with `test_rule()`.
#'
#' A response profile is the vector \eqn{t = (t_1, \dots, t_K)} of the results
#' of the \eqn{K} tests, where \eqn{t_k = 1} if test \eqn{k} is positive and
#' \eqn{t_k = 0} if it is negative. With \eqn{K} tests there are \eqn{2^K}
#' profiles. A rule assigns to each profile a value \eqn{r(t)}: \eqn{r(t) = 1}
#' if the rule declares the subject positive, \eqn{r(t) = 0} otherwise.
#'
#' A `test_rule` has two properties:
#'
#' - `name`: the label of the rule in the results;
#' - `positive`: a function of the matrix of response profiles that returns,
#'   for each profile, whether the rule declares it positive.
#'
#' The matrix of response profiles has one row per profile (\eqn{2^K} rows) and
#' one column per test, named with the test names. Its entries are `0`
#' (negative) or `1` (positive). The function in `positive` must return one
#' value per row: `TRUE` (or `1`) if the profile is positive, `FALSE` (or `0`)
#' otherwise. It may refer to tests by their names, for example
#' `profiles[, "ct"]`.
#'
#' @param name Character string, the label of the rule in the results. Use a
#'   single non-missing string.
#' @param positive Function with one argument, the matrix of response profiles
#'   described above. It must return one logical or numeric value per row of
#'   the matrix, and it is applied to the matrix of the data at the moment the
#'   rule is used.
#' @return An object of class `test_rule`.
#' @seealso [rule_accuracy()], [compare_tests()]
#' @examples
#' # Positive if the tests named "ct" and "pcr" are both positive
#' both_ct_pcr <- test_rule(
#'   name = "ct and pcr positive",
#'   positive = function(profiles) profiles[, "ct"] == 1 & profiles[, "pcr"] == 1
#' )
#' both_ct_pcr@name
#' @export
test_rule <- S7::new_class(
  "test_rule",
  properties = list(name = S7::class_character, positive = S7::class_function)
)

#' Single test rule
#'
#' A rule that declares a subject positive when one given test is positive,
#' whatever the results of the other tests.
#'
#' A response profile is the vector \eqn{t = (t_1, \dots, t_K)} of the results
#' of the \eqn{K} tests, where \eqn{t_k = 1} if test \eqn{k} is positive and
#' \eqn{t_k = 0} if it is negative. With \eqn{K} tests there are \eqn{2^K}
#' profiles.
#'
#' `single_test(test)` gives the rule
#' \deqn{r(t) = t_k,}
#' where \eqn{k} is the test named `test`. The rule declares positive the
#' profiles in which that test is positive, and ignores the results of the
#' other tests. Its label in the results is the name of the test.
#'
#' Functions that accept rules also accept a single test name in place of a
#' rule: the name is then used as a rule of this kind.
#'
#' The test name is checked when the rule is used, not when it is created: if
#' the data do not contain a test with that name, an error is raised.
#'
#' @param test Character string, the name of one test. It must be a single
#'   non-missing string that matches a test name in the data. There is no
#'   default.
#' @return An object of class `test_rule`.
#' @seealso [rule_accuracy()], [compare_tests()]
#' @examples
#' single_test("ct")
#' @export
single_test <- function(test) {
  if (!is.character(test) || length(test) != 1L || is.na(test)) {
    stop("`test` must be the name of a test.", call. = FALSE)
  }
  test_rule(name = test, positive = function(profiles) {
    if (!test %in% colnames(profiles)) {
      stop("Unknown test `", test, "`.", call. = FALSE)
    }
    profiles[, test] == 1L
  })
}

#' Any positive rule (OR rule)
#'
#' A rule that declares a subject positive when at least one test is positive.
#'
#' A response profile is the vector \eqn{t = (t_1, \dots, t_K)} of the results
#' of the \eqn{K} tests, where \eqn{t_k = 1} if test \eqn{k} is positive and
#' \eqn{t_k = 0} if it is negative. With \eqn{K} tests there are \eqn{2^K}
#' profiles.
#'
#' `any_positive()` gives the rule
#' \deqn{r(t) = 1 \textrm{ if } \sum_{k=1}^{K} t_k \ge 1, \quad r(t) = 0 \textrm{ otherwise.}}
#' The only profile declared negative is the one in which all the tests are
#' negative, \eqn{t = (0, \dots, 0)}. The label of the rule in the results is
#' `"any positive"`.
#'
#' This function has no arguments.
#'
#' @return An object of class `test_rule`.
#' @seealso [rule_accuracy()], [compare_tests()]
#' @examples
#' any_positive()
#' @export
any_positive <- function() {
  test_rule(name = "any positive", positive = function(profiles) rowSums(profiles) >= 1L)
}

#' All positive rule (AND rule)
#'
#' A rule that declares a subject positive only when all the tests are
#' positive.
#'
#' A response profile is the vector \eqn{t = (t_1, \dots, t_K)} of the results
#' of the \eqn{K} tests, where \eqn{t_k = 1} if test \eqn{k} is positive and
#' \eqn{t_k = 0} if it is negative. With \eqn{K} tests there are \eqn{2^K}
#' profiles.
#'
#' `all_positive()` gives the rule
#' \deqn{r(t) = 1 \textrm{ if } \sum_{k=1}^{K} t_k = K, \quad r(t) = 0 \textrm{ otherwise.}}
#' The only profile declared positive is the one in which all the tests are
#' positive, \eqn{t = (1, \dots, 1)}. The label of the rule in the results is
#' `"all positive"`.
#'
#' This function has no arguments.
#'
#' @return An object of class `test_rule`.
#' @seealso [rule_accuracy()], [compare_tests()]
#' @examples
#' all_positive()
#' @export
all_positive <- function() {
  test_rule(name = "all positive",
            positive = function(profiles) rowSums(profiles) == ncol(profiles))
}

#' Minimum positive rule
#'
#' A rule that declares a subject positive when at least `k` of the tests are
#' positive.
#'
#' A response profile is the vector \eqn{t = (t_1, \dots, t_K)} of the results
#' of the \eqn{K} tests, where \eqn{t_k = 1} if test \eqn{k} is positive and
#' \eqn{t_k = 0} if it is negative. With \eqn{K} tests there are \eqn{2^K}
#' profiles.
#'
#' `min_positive(k)` gives the rule
#' \deqn{r(t) = 1 \textrm{ if } \sum_{k=1}^{K} t_k \ge k, \quad r(t) = 0 \textrm{ otherwise.}}
#' With `k = 1` the rule declares positive every profile with at least one
#' positive test; with `k` equal to the number of tests it declares positive
#' only the profile with all the tests positive. The label of the rule in the
#' results is `"at least k positive"`, with the value of `k`.
#'
#' @param k Single whole number, at least 1, given as a number (for example
#'   `2`). It must not exceed the number of tests in the data: this is checked
#'   when the rule is used, and a larger value raises an error. There is no
#'   default.
#' @return An object of class `test_rule`.
#' @seealso [rule_accuracy()], [compare_tests()]
#' @examples
#' min_positive(2)
#' @export
min_positive <- function(k) {
  if (!is_scalar_number(k) || k < 1 || k != round(k)) {
    stop("`k` must be a positive whole number.", call. = FALSE)
  }
  test_rule(name = paste("at least", k, "positive"), positive = function(profiles) {
    if (k > ncol(profiles)) {
      stop("min_positive(", k, ") needs at least ", k, " tests.", call. = FALSE)
    }
    rowSums(profiles) >= k
  })
}

S7::method(print_generic, test_rule) <- function(x, ...) {
  cat("<test_rule> ", x@name, "\n", sep = "")
  invisible(x)
}

# Rules (a rule, a test name, or a list of them) as a list of test_rule objects.
as_rules <- function(rules, argument = "rules") {
  if (is.character(rules)) rules <- lapply(rules, single_test)
  if (S7::S7_inherits(rules, test_rule)) rules <- list(rules)
  if (!is.list(rules) || length(rules) == 0L) {
    stop("`", argument, "` must be a rule, the name of a test, or a list of them.",
         call. = FALSE)
  }
  lapply(rules, function(r) {
    if (is.character(r) && length(r) == 1L) r <- single_test(r)
    if (!S7::S7_inherits(r, test_rule)) {
      stop("`", argument, "` must contain rules, such as any_positive(), or test names.",
           call. = FALSE)
    }
    r
  })
}

# Profiles x rules 0/1 matrix: entry (j, r) is 1 if rule r declares profile j positive.
rule_matrix <- function(rules, data) {
  m <- vapply(rules, function(r) as.numeric(r@positive(data@profiles)),
              numeric(nrow(data@profiles)))
  m <- matrix(m, nrow = nrow(data@profiles))
  colnames(m) <- vapply(rules, function(r) r@name, character(1))
  m
}
