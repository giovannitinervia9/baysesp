test_that("gauss_hermite() integrates against the standard normal", {
  gh <- gauss_hermite(32L)
  expect_equal(sum(gh$weights), 1)
  expect_equal(sum(gh$weights * gh$nodes), 0, tolerance = 1e-10)
  expect_equal(sum(gh$weights * gh$nodes^2), 1)
  expect_equal(sum(gh$weights * gh$nodes^4), 3)
  # the closed form of the sensitivity: E Phi(a + bZ) = Phi(a / sqrt(1 + b^2))
  expect_equal(sum(gh$weights * stats::pnorm(0.8 + 1.5 * gh$nodes)),
               stats::pnorm(0.8 / sqrt(1 + 1.5^2)), tolerance = 1e-6)
})

test_that("random_effects() has the documented defaults and checks its priors", {
  m <- random_effects()
  expect_true(S7::S7_inherits(m, accuracy_model))
  expect_identical(c(m@sensitivity_prior@shape1, m@sensitivity_prior@shape2), c(1, 1))
  expect_error(random_effects(specificity_prior = prior_normal(0, 1)), "specificity_prior")
  expect_identical(m@loadings_prior@sd, 1)
  expect_output(print(m), "random effects model")
  expect_error(random_effects(loadings_prior = prior_beta()), "loadings_prior")
  expect_true(m@positive_loadings)
  expect_output(print(random_effects(positive_loadings = FALSE)), "first loading positive")
  expect_error(random_effects(positive_loadings = NA), "positive_loadings")
})

test_that("stan_data() translates the random effects model", {
  d <- example_fit(verified = c(FALSE, rep(TRUE, 7)))@data
  s <- stan_data(random_effects(), mar(), d)
  expect_identical(c(s$J, s$K, s$Q), c(8L, 3L, 32L))
  expect_identical(s$profiles[6, ], c(1, 0, 1))
  expect_equal(sum(exp(s$log_weights)), 1)
  expect_identical(c(s$mnar, s$n_beta), c(0L, 0L))
  expect_identical(s$accuracy_family, matrix(4L, 2, 3))   # beta on the probit
  expect_identical(s$positive_loadings, 1L)
  expect_identical(stan_data(random_effects(positive_loadings = FALSE), mar(), d)$positive_loadings, 0L)
  m <- random_effects(verification = mnar())
  expect_warning(s <- stan_data(m, m@verification, d), "no verified subject")
  expect_identical(c(s$mnar, s$n_beta), c(1L, 1L))
  expect_identical(stan_program(m), "random_effects")

  x <- data.frame(T = c(0, 1), V = 1, D = c(0, 1))
  expect_error(stan_data(random_effects(), mar(), baysesp_data(x, "T", "V", "D")),
               "at least two tests")
})

test_that("the random effects model fits conditionally independent data", {
  skip_on_cran()
  set.seed(21)
  s <- simulate_study(n = 10000, prevalence = 0.2, sensitivity = c(0.85, 0.75, 0.8),
                      specificity = c(0.9, 0.95, 0.85),
                      verification_prob = c(0.2, rep(0.8, 6), 1))
  d <- baysesp_data(s, c("T1", "T2", "T3"), "verified", "disease")
  fit <- baysesp(d, random_effects(), seed = 22)
  theta <- core_draws(fit)
  expect_equal(rowSums(theta$theta_diseased), rep(1, 4000), tolerance = 1e-8)
  # same data, conditional independence: the two models agree
  a <- test_accuracy(fit)
  b <- test_accuracy(baysesp(d, loglinear(), seed = 23))
  expect_equal(a$mean, b$mean, tolerance = 0.02)
  expect_output(print(fit), "random_effects model, MAR verification")

  # MNAR runs and keeps the verification parameter
  fit_mnar <- baysesp(d, random_effects(verification = mnar()), seed = 24)
  expect_true("beta[1]" %in% posterior::variables(fit_mnar@draws))
})

test_that("the random effects model recovers the accuracy under dependence", {
  skip_on_cran()
  # data from the random effects model itself, with dependent tests
  set.seed(31)
  n <- 10000
  a1 <- c(1.5, 1, 1.2)
  b1 <- c(1, 0.8, 1.2)
  a0 <- c(-1.5, -2, -1.2)
  b0 <- c(0.6, 0.5, 0.8)
  disease <- stats::rbinom(n, 1, 0.2)
  z <- stats::rnorm(n)
  a <- ifelse(disease == 1, 1, 0) %o% a1 + ifelse(disease == 0, 1, 0) %o% a0
  b <- ifelse(disease == 1, 1, 0) %o% b1 + ifelse(disease == 0, 1, 0) %o% b0
  tests <- matrix(stats::rbinom(length(a), 1, stats::pnorm(a + b * z)), n)
  profile <- 1 + tests %*% 2^(0:2)
  verified <- stats::rbinom(n, 1, c(0.2, rep(0.8, 6), 1)[profile])
  s <- data.frame(tests, verified = verified, disease = ifelse(verified == 1, disease, NA))
  d <- baysesp_data(s, c("X1", "X2", "X3"), "verified", "disease")
  fit <- baysesp(d, random_effects(), seed = 32)
  expect_identical(rstan::get_num_divergent(fit@stan_fit), 0L)
  expect_lt(max(posterior::summarise_draws(fit@draws, "rhat")$rhat), 1.01)

  sick <- disease == 1
  truth <- c(prevalence = mean(sick),
             stats::setNames(colMeans(tests[sick, ]), paste0("Se[X", 1:3, "]")),
             stats::setNames(colMeans(1 - tests[!sick, ]), paste0("Sp[X", 1:3, "]")))
  acc <- test_accuracy(fit, level = 0.99)
  for (q in names(truth)) {
    i <- acc$quantity == q
    expect_true(acc$lower[i] <= truth[[q]] && truth[[q]] <= acc$upper[i], label = q)
  }
  # the dependence matters: conditional independence misses the AND rule
  and_true <- mean(rowSums(tests[sick, ]) == 3)
  re_and <- rule_accuracy(fit, all_positive(), level = 0.99)
  ci_and <- rule_accuracy(baysesp(d, loglinear(), seed = 33), all_positive(), level = 0.99)
  i <- re_and$quantity == "Se[all positive]"
  expect_true(re_and$lower[i] <= and_true && and_true <= re_and$upper[i])
  expect_false(ci_and$lower[i] <= and_true && and_true <= ci_and$upper[i])
})

test_that("with only the first loading positive the model still runs", {
  skip_on_cran()
  set.seed(21)
  s <- simulate_study(n = 2000, prevalence = 0.2, sensitivity = c(0.85, 0.75, 0.8),
                      specificity = c(0.9, 0.95, 0.85),
                      verification_prob = c(0.2, rep(0.8, 6), 1))
  d <- baysesp_data(s, c("T1", "T2", "T3"), "verified", "disease")
  fit <- suppressWarnings(baysesp(d, random_effects(positive_loadings = FALSE), seed = 25))
  expect_identical(posterior::ndraws(fit@draws), 4000L)
})
