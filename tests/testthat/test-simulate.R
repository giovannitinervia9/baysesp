test_that("simulate_study returns individual data", {
  set.seed(1)
  d <- simulate_study(n = 50, prevalence = 0.3, sensitivity = c(A = 0.9, B = 0.8),
                      specificity = c(0.9, 0.9), verification_prob = 0.5)
  expect_identical(names(d), c("A", "B", "verified", "disease"))
  expect_identical(nrow(d), 50L)
  expect_true(all(is.na(d$disease[d$verified == 0])))
  expect_false(anyNA(d$disease[d$verified == 1]))
  expect_identical(names(simulate_study(5, 0.1, 0.9, 0.9, 1)), c("T1", "verified", "disease"))
})

test_that("simulate_study reproduces the parameters with many subjects", {
  set.seed(2)
  n <- 2e5
  d <- simulate_study(n = n, prevalence = 0.2, sensitivity = c(0.9, 0.6),
                      specificity = c(0.8, 0.95), verification_prob = 1)
  tol <- 0.01
  expect_equal(mean(d$disease), 0.2, tolerance = tol)
  sick <- d$disease == 1
  expect_equal(colMeans(d[sick, c("T1", "T2")]), c(T1 = 0.9, T2 = 0.6), tolerance = tol)
  expect_equal(colMeans(1 - d[!sick, c("T1", "T2")]), c(T1 = 0.8, T2 = 0.95),
               tolerance = tol)
  # conditional independence: P(T1 = 1, T2 = 1 | D = 1) = 0.9 * 0.6
  expect_equal(mean(d$T1[sick] * d$T2[sick]), 0.54, tolerance = tol)
})

test_that("verification follows the profile and the odds ratio", {
  set.seed(3)
  # only profiles (1, 0) and (1, 1) are verified
  d <- simulate_study(n = 2000, prevalence = 0.3, sensitivity = c(0.8, 0.7),
                      specificity = c(0.8, 0.8), verification_prob = c(0, 1, 0, 1))
  expect_identical(d$verified, as.integer(d$T1 == 1))
  b <- baysesp_data(d, c("T1", "T2"), "verified", "disease")
  expect_identical(b@never_verified, c(TRUE, FALSE, TRUE, FALSE))

  # MNAR: odds of verification 3 times higher in the diseased
  set.seed(4)
  n <- 2e5
  d <- simulate_study(n = n, prevalence = 0.5, sensitivity = 0.5, specificity = 0.5,
                      verification_prob = 0.2, odds_ratio = 3)
  verified_sick <- sum(d$disease == 1, na.rm = TRUE) / (n / 2)
  expect_equal(verified_sick, stats::plogis(stats::qlogis(0.2) + log(3)), tolerance = 0.01)
})

test_that("simulate_study rejects invalid arguments", {
  expect_error(simulate_study(0, 0.1, 0.9, 0.9, 1), "`n`")
  expect_error(simulate_study(10, 1.1, 0.9, 0.9, 1), "`prevalence`")
  expect_error(simulate_study(10, 0.1, c(0.9, 0.8), 0.9, 1), "same length")
  expect_error(simulate_study(10, 0.1, 0.9, 0.9, c(1, 1, 1)), "one per profile")
  expect_error(simulate_study(10, 0.1, 0.9, 0.9, 1, odds_ratio = 0), "`odds_ratio`")
  expect_error(simulate_study(10, 0.1, c(disease = 0.9), 0.9, 1), "cannot be called")
})
