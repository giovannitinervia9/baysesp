test_that("mnar() has the documented defaults and checks its priors", {
  v <- mnar()
  expect_true(S7::S7_inherits(v, verification_mechanism))
  expect_identical(v@odds_ratio_prior@sd, 1)
  expect_identical(v@intercepts_prior@sd, 1.5)
  expect_false(v@by_profile)
  expect_output(print(mnar(by_profile = TRUE)), "one per profile")
  expect_error(mnar(odds_ratio_prior = prior_beta()), "odds_ratio_prior")
  expect_error(mnar(by_profile = NA), "by_profile")
})

test_that("stan_data() translates mnar() in the saturated model", {
  x <- data.frame(T1 = c(0, 1, 0, 1), V = c(1, 1, 0, 0), D = c(0, 1, NA, NA))
  d <- baysesp_data(x, "T1", "V", "D")
  m <- saturated(verification = mnar(prior_student(3, 0.5, 1)))
  s <- stan_data(m, m@verification, d)
  expect_identical(s$n_beta, 1L)
  expect_identical(s$beta_family, 2L)
  expect_equal(s$beta_params, c(3, 0.5, 1))
  m <- saturated(verification = mnar(by_profile = TRUE))
  expect_identical(stan_data(m, m@verification, d)$n_beta, 2L)
  expect_identical(stan_data(saturated(), mar(), d)$n_beta, 0L)
})

test_that("with odds ratio 1 and the same priors the analysis reproduces the fit", {
  fit <- example_fit()
  s <- sensitivity_analysis(fit, odds_ratios = 1)
  acc <- test_accuracy(fit)
  expect_identical(names(s), c("odds_ratio", "unverified_prior", "quantity", "mean",
                               "median", "lower", "upper"))
  expect_equal(s$mean, acc$mean)
  expect_equal(sensitivity_analysis(fit)$mean, acc$mean)
  expect_true(all(s$odds_ratio == 1))
  expect_true(all(is.na(s$unverified_prior)))
})

test_that("the odds ratio changes the proportion of diseased among the unverified", {
  fit <- example_fit()
  s <- sensitivity_analysis(fit, odds_ratios = c(1, 3), rules = any_positive())
  expect_identical(unique(s$odds_ratio), c(1, 3))
  # by hand: psi = v phi1 + (1 - v) expit(logit(phi1) - log(3)), v = 0.5
  omega <- 0.1 * example_theta1 + 0.9 * example_theta0
  phi1 <- 0.1 * example_theta1 / omega
  psi <- 0.5 * phi1 + 0.5 * stats::plogis(stats::qlogis(phi1) - log(3))
  tpr_or <- sum((omega * psi)[-1]) / sum(omega * psi)
  expect_equal(s$mean[s$odds_ratio == 3 & s$quantity == "Se[any positive]"], tpr_or)
  expect_error(sensitivity_analysis(fit, odds_ratios = 0), "`odds_ratios`")
})

test_that("new priors on the never-verified profiles are sampled", {
  fit <- example_fit(verified = c(FALSE, rep(TRUE, 7)))
  set.seed(1)
  s <- sensitivity_analysis(fit, unverified_priors = list(prior_beta(0.01, 1e6),
                                                          prior_beta(0.5, 1e6)))
  expect_identical(length(unique(s$unverified_prior)), 2L)
  expect_match(unique(s$unverified_prior), "^Beta[(]mean = 0[.](01|5), precision")
  omega <- 0.1 * example_theta1 + 0.9 * example_theta0
  c0 <- sum(0.1 * example_theta1[-1])
  prev <- s$mean[s$quantity == "prevalence"]
  expect_equal(prev, c0 + omega[1] * c(0.01, 0.5), tolerance = 1e-3)

  expect_error(sensitivity_analysis(example_fit(), unverified_priors = prior_beta()),
               "every profile has verified subjects")
  expect_error(sensitivity_analysis(fit, unverified_priors = list(prior_dirichlet(1))),
               "dirichlet prior")
})

test_that("identification functions take a fixed odds ratio", {
  fit <- example_fit(verified = c(FALSE, rep(TRUE, 7)))
  b1 <- accuracy_bounds(fit)
  expect_equal(accuracy_bounds(fit, odds_ratio = 1), b1)
  b3 <- accuracy_bounds(fit, odds_ratio = 3)
  # with odds ratio 3 the unverified are less often diseased: lower prevalence bounds
  expect_true(all(b3$lower_bound[1] < b1$lower_bound[1]))
  omega <- 0.1 * example_theta1 + 0.9 * example_theta0
  phi1 <- 0.1 * example_theta1 / omega
  psi <- 0.5 * phi1 + 0.5 * stats::plogis(stats::qlogis(phi1) - log(3))
  expect_equal(b3$lower_bound[1], sum((omega * psi)[-1]))
  expect_equal(b3$upper_bound[1], sum((omega * psi)[-1]) + omega[1])
  expect_identical(nrow(accuracy_ratios(fit, odds_ratio = 3)), 6L)

  mnar_fit <- fit
  mnar_fit@model <- saturated(verification = mnar())
  expect_error(accuracy_bounds(mnar_fit), "give it with `odds_ratio`")
  expect_no_error(accuracy_bounds(mnar_fit, odds_ratio = 2))
  expect_error(accuracy_bounds(fit, odds_ratio = -1), "`odds_ratio`")
})
