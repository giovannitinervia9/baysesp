test_that("prior_beta stores mean and precision and derives the shapes", {
  p <- prior_beta(mean = 0.2, precision = 10)
  expect_true(S7::S7_inherits(p, prior_distribution))
  expect_identical(p@family, "beta")
  expect_identical(p@support, "probability")
  expect_equal(p@shape1, 2)
  expect_equal(p@shape2, 8)
  u <- prior_beta()
  expect_equal(c(u@shape1, u@shape2), c(1, 1))
})

test_that("prior_beta rejects invalid values", {
  expect_error(prior_beta(mean = 0), "@mean")
  expect_error(prior_beta(mean = 1.2), "@mean")
  expect_error(prior_beta(mean = c(0.2, 0.3)), "@mean")
  expect_error(prior_beta(mean = 0.5, precision = -1), "@precision")
})

test_that("the shapes of prior_beta are read-only", {
  p <- prior_beta(0.3, 5)
  expect_error(p@shape1 <- 2)
})

test_that("elicit_beta matches the requested upper quantile", {
  p <- elicit_beta(mean = 0.005, upper = 0.015)
  expect_true(S7::S7_inherits(p, prior_beta))
  expect_equal(p@mean, 0.005)
  expect_equal(stats::pbeta(0.015, p@shape1, p@shape2), 0.95, tolerance = 1e-6)
  # the regular Beta, not the U-shaped one with mass at the extremes
  expect_gt(p@shape1, 0.5)

  q <- elicit_beta(mean = 0.8, upper = 0.95, prob = 0.9)
  expect_equal(stats::pbeta(0.95, q@shape1, q@shape2), 0.9, tolerance = 1e-6)
})

test_that("elicit_beta rejects impossible requests", {
  expect_error(elicit_beta(mean = 0.5, upper = 0.4), "`upper`")
  expect_error(elicit_beta(mean = 1.5, upper = 0.4), "`mean`")
  expect_error(elicit_beta(mean = 0.5, upper = 0.6, prob = 0.3), "`prob`")
})

test_that("prior_normal, prior_student and prior_dirichlet validate their parameters", {
  expect_identical(prior_normal()@support, "real")
  expect_error(prior_normal(0, 0), "@sd")
  s <- prior_student(df = 3, location = 1, scale = 2)
  expect_identical(c(s@df, s@location, s@scale), c(3, 1, 2))
  expect_error(prior_student(df = -1), "@df")
  expect_error(prior_student(scale = 0), "@scale")
  d <- prior_dirichlet(c(2, 1, 1))
  expect_identical(d@support, "simplex")
  expect_error(prior_dirichlet(c(1, 0)), "@concentration")
})

test_that("the abstract class cannot be instantiated", {
  expect_error(prior_distribution(family = "x", support = "real"))
})

test_that("priors print a readable description", {
  expect_output(print(prior_beta(0.5, 2)),
                "Beta\\(mean = 0.5, precision = 2\\)\n")
  expect_output(print(prior_normal(0, 1.5)), "<prior_normal> Normal\\(mean = 0, sd = 1.5\\)")
  expect_output(print(prior_student(3, 0, 1)), "logit")
  expect_output(print(prior_dirichlet(1)), "vector of probabilities")
})

test_that("prior_code gives the Stan family code and three parameters", {
  expect_identical(prior_code(prior_normal(1, 2)), list(family = 1L, params = c(1, 2, 0)))
  expect_identical(prior_code(prior_student(3, 0, 1)), list(family = 2L, params = c(3, 0, 1)))
  expect_identical(prior_code(prior_beta(0.25, 4)), list(family = 3L, params = c(1, 3, 0)))
  expect_error(prior_code(prior_dirichlet(1)), "concentration")
})

test_that("check_prior allows only compatible scales", {
  expect_silent(check_prior(prior_beta(), "probability", "x"))
  expect_silent(check_prior(prior_normal(), "probability", "x"))
  expect_silent(check_prior(prior_dirichlet(1), "simplex", "x"))
  expect_error(check_prior(prior_beta(), "real", "x"), "beta prior")
  expect_error(check_prior(prior_normal(), "simplex", "x"), "normal prior")
  expect_error(check_prior(1, "real", "x"), "must be a prior")
})
