#' Diagnostic rules
#'
#' A rule declares a subject positive or negative from the results of the
#' tests. Its accuracy depends on the joint distribution of the tests, not only
#' on their sensitivities and specificities.
#'
#' - `single_test(test)`: positive if the test is positive;
#' - `any_positive()`: positive if at least one test is positive (OR rule);
#' - `all_positive()`: positive if all tests are positive (AND rule);
#' - `min_positive(k)`: positive if at least `k` tests are positive.
#'
#' Wherever a rule is expected, the name of a test can be given instead.
#'
#' @param test Name of a test.
#' @param k Minimum number of positive tests.
#' @param name Label of the rule in the results.
#' @param positive Function of the matrix of response profiles (one row per
#'   profile, one named column per test) that returns, for each profile,
#'   whether the rule declares it positive.
#' @return An object of class `test_rule`.
#' @seealso [rule_accuracy()], [compare_tests()]
#' @examples
#' any_positive()
#' min_positive(2)
#' @export
test_rule <- S7::new_class(
  "test_rule",
  properties = list(name = S7::class_character, positive = S7::class_function)
)

#' @rdname test_rule
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

#' @rdname test_rule
#' @export
any_positive <- function() {
  test_rule(name = "any positive", positive = function(profiles) rowSums(profiles) >= 1L)
}

#' @rdname test_rule
#' @export
all_positive <- function() {
  test_rule(name = "all positive",
            positive = function(profiles) rowSums(profiles) == ncol(profiles))
}

#' @rdname test_rule
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

S7::method(print, test_rule) <- function(x, ...) {
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
