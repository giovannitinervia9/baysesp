test_that("design_matrix() gives the matrix of Example 5.2", {
  profiles <- example_fit()@data@profiles
  X <- design_matrix(profiles, 3)
  expect_identical(colnames(X), c("A", "B", "C", "A:B", "A:C", "B:C", "A:B:C"))
  expect_identical(attr(X, "order"), c(1L, 1L, 1L, 2L, 2L, 2L, 3L))
  expect_identical(unname(X[4, ]), c(1, 1, 0, 1, 0, 0, 0))
  expect_identical(unname(X[8, ]), rep(1, 7))
  expect_identical(unname(X[6, ]), c(1, 0, 1, 0, 1, 0, 0))
  # lower orders keep the first columns
  expect_identical(unname(design_matrix(profiles, 1)[, 1:3]), unname(X[, 1:3]))
  expect_identical(ncol(design_matrix(profiles, 2)), 6L)
})

test_that("loglinear() has the documented defaults and checks its priors", {
  m <- loglinear()
  expect_true(S7::S7_inherits(m, accuracy_model))
  expect_identical(m@interaction_order, 1L)
  expect_identical(m@main_effects_prior@sd, 1.5)
  expect_identical(m@interactions_prior@sd, 1)
  expect_output(print(loglinear(2)), "log-linear model of order 2")
  expect_silent(loglinear(prevalence_prior = prior_normal(-2, 1)))
  expect_error(loglinear(main_effects_prior = prior_beta()), "main_effects_prior")
  expect_error(loglinear(interaction_order = 0), "interaction_order")
  expect_error(loglinear(interaction_order = 1.5), "interaction_order")
})

test_that("stan_data() translates the log-linear models", {
  d <- example_fit(verified = c(FALSE, rep(TRUE, 7)))@data
  s <- stan_data(loglinear(2), mar(), d)
  expect_identical(c(s$J, s$P), c(8L, 6L))
  expect_identical(s$order, as.array(c(1L, 1L, 1L, 2L, 2L, 2L)))
  expect_identical(c(s$mnar, s$n_beta), c(0L, 0L))
  expect_identical(s$prevalence_family, 3L)
  m <- loglinear(verification = mnar(by_profile = TRUE))
  s <- stan_data(m, m@verification, d)
  expect_identical(c(s$mnar, s$n_beta), c(1L, 7L))
  expect_identical(s$alpha_params, c(0, 2.5, 0))
  expect_identical(stan_program(m), "loglinear")
  expect_error(stan_data(loglinear(4), mar(), d), "cannot exceed the number of tests \\(3\\)")
})

test_that("functions for the saturated model refuse the structured models", {
  fit <- example_fit(verified = c(FALSE, rep(TRUE, 7)))
  fit@model <- loglinear()
  expect_no_error(test_accuracy(fit))
  expect_error(identification_bounds(fit), "only for the saturated model")
  expect_error(sensitivity_analysis(fit), "`odds_ratios`")
  expect_error(sensitivity_analysis(fit, odds_ratios = 2, unverified_priors = prior_beta()),
               "applies only to saturated")
})

test_that("conditional independence recovers the accuracy with a never-verified profile", {
  skip_on_cran()
  # profile (0,0) never verified; the tests are conditionally independent
  set.seed(11)
  s <- simulate_study(n = 20000, prevalence = 0.1, sensitivity = c(0.85, 0.75),
                      specificity = c(0.9, 0.95), verification_prob = c(0, 0.8, 0.8, 1))
  d <- baysesp_data(s, c("T1", "T2"), "verified", "disease")
  fit <- baysesp(d, loglinear(), seed = 12)
  expect_setequal(posterior::variables(fit@draws),
                  c("prevalence", paste0("theta_healthy[", 1:4, "]"),
                    paste0("theta_diseased[", 1:4, "]")))
  # the values in this sample: same seed, everyone verified (verification is
  # drawn last, so disease and test results are the same)
  set.seed(11)
  full <- simulate_study(n = 20000, prevalence = 0.1, sensitivity = c(0.85, 0.75),
                         specificity = c(0.9, 0.95), verification_prob = 1)
  sick <- full$disease == 1
  truth <- c(prevalence = mean(sick),
             `Se[T1]` = mean(full$T1[sick]), `Se[T2]` = mean(full$T2[sick]),
             `Sp[T1]` = mean(1 - full$T1[!sick]), `Sp[T2]` = mean(1 - full$T2[!sick]))
  acc <- test_accuracy(fit, level = 0.99)
  for (q in names(truth)) {
    i <- acc$quantity == q
    expect_true(acc$lower[i] <= truth[[q]] && truth[[q]] <= acc$upper[i], label = q)
  }
  # identified by the structure, not by the data
  expect_false(any(acc$identified[!startsWith(acc$quantity, "PPV")]))
  expect_output(print(fit), "loglinear model, MAR verification")
})

test_that("the log-linear model of order K agrees with the saturated model", {
  skip_on_cran()
  set.seed(13)
  s <- simulate_study(n = 20000, prevalence = 0.2, sensitivity = c(0.85, 0.75),
                      specificity = c(0.9, 0.95), verification_prob = c(0.3, 0.8, 0.8, 1))
  d <- baysesp_data(s, c("T1", "T2"), "verified", "disease")
  a <- test_accuracy(baysesp(d, loglinear(2), seed = 14))
  b <- test_accuracy(baysesp(d, saturated(), seed = 15))
  expect_equal(a$mean, b$mean, tolerance = 0.01)
})

test_that("MNAR in the log-linear model: fit and refits at fixed odds ratios", {
  skip_on_cran()
  set.seed(16)
  s <- simulate_study(n = 5000, prevalence = 0.15, sensitivity = c(0.8, 0.7),
                      specificity = c(0.9, 0.95), verification_prob = c(0.2, 0.7, 0.7, 0.9),
                      odds_ratio = exp(1))
  d <- baysesp_data(s, c("T1", "T2"), "verified", "disease")
  fit <- baysesp(d, loglinear(verification = mnar()), seed = 17)
  expect_true("beta[1]" %in% posterior::variables(fit@draws))
  expect_output(print(fit), "MNAR verification")

  sens <- sensitivity_analysis(fit, odds_ratios = c(1, exp(1)), seed = 18)
  expect_identical(unique(sens$odds_ratio), c(1, exp(1)))
  # odds ratio 1 is the MAR model
  mar_acc <- test_accuracy(baysesp(d, loglinear(), seed = 19))
  expect_equal(sens$mean[sens$odds_ratio == 1], mar_acc$mean, tolerance = 0.02)
})
