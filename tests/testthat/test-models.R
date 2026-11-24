test_that("saturated() has the documented defaults", {
  m <- saturated()
  expect_true(S7::S7_inherits(m, accuracy_model))
  expect_true(S7::S7_inherits(m@verification, mar))
  expect_identical(m@frequency_prior@concentration, 1)
  expect_equal(c(m@verified_prior@shape1, m@verified_prior@shape2), c(1, 1))
  expect_equal(c(m@profile_prevalence_prior@shape1, m@profile_prevalence_prior@shape2),
               c(1, 1))
  expect_null(m@unverified_prevalence_prior)
  expect_output(print(m), "unverified_prevalence_prior: none")
  expect_error(accuracy_model(), "abstract")
  expect_error(verification_mechanism(), "abstract")
})

test_that("saturated() accepts only compatible priors", {
  expect_silent(saturated(verified_prior = prior_normal(0, 1.5),
                          unverified_prevalence_prior = list(prior_beta(0.1, 20),
                                                             prior_student())))
  expect_error(saturated(verified_prior = prior_dirichlet(1)), "verified_prior")
  expect_error(saturated(frequency_prior = prior_beta()), "frequency_prior")
  expect_error(saturated(unverified_prevalence_prior = list(prior_beta(), 1)),
               "unverified_prevalence_prior\\[\\[2\\]\\]")
  expect_error(saturated(verification = "mar"), "verification")
})

test_that("stan_data() translates the saturated model with MAR", {
  x <- data.frame(T1 = c(0, 1, 0, 1, 1), T2 = c(0, 0, 1, 1, 1),
                  V = c(0, 1, 1, 1, 0), D = c(NA, 1, 0, 1, NA))
  d <- baysesp_data(x, c("T1", "T2"), "V", "D")
  m <- saturated(unverified_prevalence_prior = prior_beta(0.1, 10))
  s <- stan_data(m, m@verification, d)
  expect_identical(s$J, 4L)
  expect_identical(s$n1, c(0L, 1L, 0L, 1L))
  expect_identical(s$n0, c(0L, 0L, 1L, 0L))
  expect_identical(s$u, c(1L, 0L, 0L, 1L))
  expect_identical(s$frequency_concentration, rep(1, 4))
  expect_identical(s$verified_family, 3L)
  expect_equal(s$verified_params, c(1, 1, 0))
  expect_identical(s$J0, 1L)
  expect_identical(s$q_family, as.array(3L))
  expect_equal(s$q_params, matrix(c(1, 9, 0), 1, 3))
  expect_identical(s$mnar, 0L)
  expect_identical(stan_program(m), "saturated")

  # without never-verified profiles the q arrays are empty
  d_all <- baysesp_data(transform(x, V = 1, D = c(0, 1, 0, 1, 1)), c("T1", "T2"), "V", "D")
  s_all <- stan_data(m, m@verification, d_all)
  expect_identical(s_all$J0, 0L)
  expect_identical(dim(s_all$q_params), c(0L, 3L))
})

test_that("stan_data() requires a prior for each never-verified profile", {
  x <- data.frame(T1 = c(0, 1, 0, 1), T2 = c(0, 0, 1, 1),
                  V = c(0, 1, 0, 1), D = c(NA, 1, NA, 0))
  d <- baysesp_data(x, c("T1", "T2"), "V", "D")
  m <- saturated()
  expect_error(stan_data(m, m@verification, d), "profiles \\(0,0\\) \\(0,1\\)")
  m <- saturated(unverified_prevalence_prior = list(prior_beta(), prior_beta(), prior_beta()))
  expect_error(stan_data(m, m@verification, d), "one prior per never-verified profile \\(2\\)")
  m <- saturated(frequency_prior = prior_dirichlet(c(1, 2)),
                 unverified_prevalence_prior = prior_beta())
  expect_error(stan_data(m, m@verification, d), "one per profile \\(4\\)")
})

test_that("test_accuracy() computes the accuracy from the profile distributions", {
  # every profile verified except (0,0,0)
  fit <- example_fit(verified = c(FALSE, rep(TRUE, 7)))

  acc <- test_accuracy(fit)
  get <- function(q) acc$mean[acc$quantity == q]
  expect_identical(acc$quantity[1:5], c("prevalence", "Se[A]", "Sp[A]", "PPV[A]", "NPV[A]"))
  expect_identical(nrow(acc), 13L)
  expect_equal(get("prevalence"), 0.1)
  expect_equal(get("Se[A]"), 0.80)
  expect_equal(get("Sp[A]"), 0.90)
  expect_equal(get("PPV[A]"), 0.1 * 0.8 / (0.1 * 0.8 + 0.9 * 0.1))
  expect_equal(get("NPV[A]"), 0.9 * 0.9 / (0.9 * 0.9 + 0.1 * 0.2))
  expect_equal(get("Se[C]"), 0.03 + 0.12 + 0.08 + 0.50)
  expect_identical(acc$lower, acc$mean)
  # with (0,0,0) never verified only the PPVs are identified
  expect_identical(acc$identified, c(FALSE, rep(c(FALSE, FALSE, TRUE, FALSE), 3)))
  expect_error(test_accuracy(fit, level = 1), "`level`")
  expect_error(test_accuracy(1), "`fit`")
})
