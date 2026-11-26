# Rule vectors and accuracy of Example 4.5 of the specification (K = 3).

test_that("rules give the vectors of Example 4.5", {
  d <- example_fit()@data
  m <- rule_matrix(as_rules(list("A", any_positive(), all_positive(), min_positive(2))), d)
  expect_identical(colnames(m), c("A", "any positive", "all positive", "at least 2 positive"))
  expect_identical(m[, 1], c(0, 1, 0, 1, 0, 1, 0, 1))
  expect_identical(m[, 2], c(0, 1, 1, 1, 1, 1, 1, 1))
  expect_identical(m[, 3], c(0, 0, 0, 0, 0, 0, 0, 1))
  expect_identical(m[, 4], c(0, 0, 0, 1, 0, 1, 1, 1))
  expect_output(print(any_positive()), "<test_rule> any positive")

  expect_error(rule_matrix(as_rules("Z"), d), "Unknown test `Z`")
  expect_error(rule_matrix(as_rules(min_positive(4)), d), "at least 4 tests")
  expect_error(min_positive(0), "`k`")
  expect_error(single_test(1), "`test`")
  expect_error(as_rules(list(1)), "must contain rules")
})

test_that("rule_accuracy() reproduces the table of Example 4.5", {
  acc <- rule_accuracy(example_fit(), list("A", any_positive(), all_positive(),
                                            min_positive(2)))
  get <- function(q) acc$mean[acc$quantity == q]
  expect_identical(acc$quantity[1:4], c("Se[A]", "Sp[A]", "PPV[A]", "NPV[A]"))
  expect_equal(get("Se[A]"), 0.80)
  expect_equal(get("Sp[A]"), 0.90)
  expect_equal(get("Se[any positive]"), 0.95)
  expect_equal(get("Sp[any positive]"), 0.80)
  expect_equal(get("Se[all positive]"), 0.50)
  expect_equal(get("Sp[all positive]"), 0.99)
  expect_equal(get("Se[at least 2 positive]"), 0.80)
  expect_equal(get("Sp[at least 2 positive]"), 0.95)
  expect_equal(get("NPV[all positive]"), 0.9 * 0.99 / (0.9 * 0.99 + 0.1 * 0.5))
  expect_true(all(acc$identified))
})

test_that("compare_tests() gives differences, ratios and probabilities", {
  cmp <- compare_tests(example_fit(), "A", any_positive())
  expect_identical(cmp$quantity,
                   c("Se[A] - Se[any positive]", "Sp[A] - Sp[any positive]",
                     "Youden[A] - Youden[any positive]",
                     "TPR[A] / TPR[any positive]", "FPR[A] / FPR[any positive]"))
  expect_equal(cmp$mean, c(-0.15, 0.10, -0.05, 0.80 / 0.95, 0.5))
  expect_equal(cmp$probability, c(0, 1, 0, 0, 0))
  expect_error(compare_tests(example_fit(), c("A", "B"), "C"), "one test or rule each")
})

test_that("identification follows Proposition 7.1", {
  # profile (0,0,0) never verified: every rule declares it negative, so the PPVs
  # and the ratios are identified, nothing else is
  fit <- example_fit(verified = c(FALSE, rep(TRUE, 7)))
  acc <- rule_accuracy(fit, list("A", any_positive()))
  expect_identical(acc$identified[acc$quantity %in% c("PPV[A]", "PPV[any positive]")],
                   c(TRUE, TRUE))
  expect_false(any(acc$identified[!grepl("^PPV", acc$quantity)]))
  cmp <- compare_tests(fit, "A", "B")
  expect_identical(cmp$identified, c(FALSE, FALSE, FALSE, TRUE, TRUE))

  # profile (1,0,0) never verified: test A and the OR rule declare it positive
  fit <- example_fit(verified = c(TRUE, FALSE, rep(TRUE, 6)))
  expect_identical(compare_tests(fit, "A", "B")$identified, rep(FALSE, 5))
  expect_identical(compare_tests(fit, "B", "C")$identified, c(FALSE, FALSE, FALSE, TRUE, TRUE))
})

test_that("accuracy_bounds() gives the bounds of Proposition 7.1", {
  fit <- example_fit(verified = c(FALSE, rep(TRUE, 7)))
  b <- accuracy_bounds(fit, rules = list("A", any_positive()))
  p <- 0.1
  omega0 <- p * example_theta1[1] + (1 - p) * example_theta0[1]
  c0 <- p * (1 - example_theta1[1])
  A <- p * 0.80
  B <- (1 - p) * 0.10
  expect_identical(b$quantity, c("prevalence", "Se[A]", "Sp[A]",
                                 "Se[any positive]", "Sp[any positive]"))
  expect_equal(b$lower_bound[1:3], c(c0, A / (c0 + omega0), 1 - B / (1 - c0 - omega0)))
  expect_equal(b$upper_bound[1:3], c(c0 + omega0, A / c0, 1 - B / (1 - c0)))
  expect_equal(b$lower, b$lower_bound)
  # the true values are inside the bounds
  expect_true(b$lower_bound[2] <= 0.80 && 0.80 <= b$upper_bound[2])
  expect_true(b$lower_bound[1] <= p && p <= b$upper_bound[1])

  # a rule that declares a never-verified profile positive has no bounds
  fit <- example_fit(verified = c(TRUE, FALSE, rep(TRUE, 6)))
  b <- accuracy_bounds(fit)
  expect_true(all(is.na(b[b$quantity %in% c("Se[A]", "Sp[A]"), -1])))
  expect_false(anyNA(b[b$quantity == "Se[B]", -1]))
})

test_that("accuracy_ratios() compares the eligible rules", {
  fit <- example_fit(verified = c(TRUE, FALSE, rep(TRUE, 6)))
  q <- accuracy_ratios(fit)
  expect_identical(q$quantity, c("TPR[B] / TPR[C]", "FPR[B] / FPR[C]"))
  expect_equal(q$mean, c((0.04 + 0.10 + 0.08 + 0.50) / (0.03 + 0.12 + 0.08 + 0.50),
                         (0.05 + 0.02 + 0.01 + 0.01) / (0.04 + 0.01 + 0.01 + 0.01)))
  expect_identical(names(q), c("quantity", "mean", "median", "lower", "upper", "probability"))
  expect_error(accuracy_ratios(fit, rules = c("A", "B")), "At least two rules")
})
