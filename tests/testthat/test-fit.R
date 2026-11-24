# Fits of saturated.stan. With MAR and every profile verified, the posterior of
# the saturated model is known in closed form (specification, section 10.2):
# omega ~ Dirichlet(1 + N), phi1 ~ Beta(1 + n1, 1 + n0), and the prevalence is
# sum(omega * phi1). The parameters of the never-verified profiles do not enter
# the likelihood, so their posterior equals their prior.

# aggregated table of Example 2.1 of the specification (K = 2, N = 1000)
example_data <- function() {
  profiles <- data.frame(T1 = c(0, 1, 0, 1), T2 = c(0, 0, 1, 1))
  agg <- rbind(
    cbind(profiles, V = 1, D = 1, n = c(1, 8, 5, 30)),
    cbind(profiles, V = 1, D = 0, n = c(49, 22, 25, 15)),
    cbind(profiles, V = 0, D = NA, n = c(750, 70, 20, 5))
  )
  baysesp_data(agg, c("T1", "T2"), "V", "D", counts = "n")
}

# exact posterior draws of the prevalence and of Se[T1]
conjugate_draws <- function(d, n = 1e5) {
  N <- d@n_diseased + d@n_healthy + d@n_unverified
  g <- matrix(stats::rgamma(n * 4, shape = 1 + N), n, 4, byrow = TRUE)
  omega <- g / rowSums(g)
  phi1 <- matrix(stats::rbeta(n * 4, 1 + d@n_diseased, 1 + d@n_healthy), n, 4,
                 byrow = TRUE)
  prevalence <- rowSums(omega * phi1)
  se1 <- rowSums((omega * phi1)[, c(2, 4)]) / prevalence
  list(prevalence = prevalence, se1 = se1)
}

expect_same_posterior <- function(x, exact) {
  expect_lt(abs(mean(x) - mean(exact)), 0.1 * sd(exact))
  expect_equal(sd(x), sd(exact), tolerance = 0.1)
}

test_that("the saturated model reproduces the conjugate posterior (rstan)", {
  skip_on_cran()
  d <- example_data()
  fit <- baysesp(d, seed = 1)
  expect_true(S7::S7_inherits(fit, baysesp_fit))
  expect_identical(fit@backend, "rstan")
  expect_setequal(posterior::variables(fit@draws),
                  c("prevalence", paste0("theta_healthy[", 1:4, "]"),
                    paste0("theta_diseased[", 1:4, "]"), paste0("omega[", 1:4, "]"),
                    paste0("v[", 1:4, "]"), paste0("phi1[", 1:4, "]")))
  expect_identical(posterior::ndraws(fit@draws), 4000L)

  set.seed(1)
  exact <- conjugate_draws(d)
  acc <- test_accuracy(fit)
  expect_true(all(acc$identified))
  draws <- core_draws(fit)
  expect_same_posterior(draws$prevalence, exact$prevalence)
  expect_same_posterior(rowSums(draws$theta_diseased[, c(2, 4)]), exact$se1)
  expect_output(print(fit), "saturated model, MAR verification, 4 chains x 1000 draws \\(rstan\\)")
})

test_that("the cmdstanr backend gives the same posterior", {
  skip_on_cran()
  skip_if_not(cmdstan_available(), "cmdstanr and CmdStan are not available")
  d <- example_data()
  fit <- baysesp(d, backend = "cmdstanr", seed = 1)
  expect_identical(fit@backend, "cmdstanr")
  set.seed(2)
  exact <- conjugate_draws(d)
  expect_same_posterior(core_draws(fit)$prevalence, exact$prevalence)
})

test_that("in never-verified profiles the posterior equals the prior", {
  skip_on_cran()
  # only subjects positive to T1 are verified: profiles (0,0) and (0,1) never are
  set.seed(3)
  s <- simulate_study(n = 2000, prevalence = 0.2, sensitivity = c(0.9, 0.7),
                      specificity = c(0.85, 0.95), verification_prob = c(0, 0.5, 0, 0.5))
  d <- baysesp_data(s, c("T1", "T2"), "verified", "disease")
  expect_error(baysesp(d), "No subject is verified in the profiles \\(0,0\\) \\(0,1\\)")

  # a Beta prior on the first profile, a normal prior on the logit for the second
  m <- saturated(unverified_prevalence_prior = list(prior_beta(0.05, 40),
                                                    prior_normal(-2, 0.5)))
  fit <- baysesp(d, m, seed = 4, iter_sampling = 2000)
  q <- posterior::as_draws_matrix(posterior::subset_draws(fit@draws, "q"))
  expect_equal(mean(q[, 1]), 0.05, tolerance = 0.05)
  expect_equal(sd(q[, 1]), sqrt(0.05 * 0.95 / 41), tolerance = 0.1)
  expect_equal(mean(stats::qlogis(q[, 2])), -2, tolerance = 0.05)
  expect_equal(sd(stats::qlogis(q[, 2])), 0.5, tolerance = 0.1)

  acc <- test_accuracy(fit)
  # the PPV of T1 involves only verified profiles
  expect_identical(acc$quantity[acc$identified], "PPV[T1]")
})

test_that("with MNAR the posterior of beta equals its prior and matches the R analysis", {
  skip_on_cran()
  set.seed(5)
  s <- simulate_study(n = 5000, prevalence = 0.2, sensitivity = c(0.85, 0.7),
                      specificity = c(0.9, 0.95), verification_prob = c(0.2, 0.7, 0.7, 0.9),
                      odds_ratio = exp(1))
  d <- baysesp_data(s, c("T1", "T2"), "verified", "disease")

  # beta ~ N(0, 1) is not updated by the data
  fit <- baysesp(d, saturated(verification = mnar()), seed = 6, iter_sampling = 2000)
  beta <- as.numeric(posterior::as_draws_matrix(posterior::subset_draws(fit@draws, "beta")))
  expect_equal(mean(beta), 0, tolerance = 0.05)
  expect_equal(sd(beta), 1, tolerance = 0.05)
  expect_false(any(test_accuracy(fit)$identified))

  # an almost fixed beta = 1 in Stan agrees with the analysis at odds ratio e in R
  fixed <- baysesp(d, saturated(verification = mnar(prior_normal(1, 0.001))), seed = 7)
  mar_fit <- baysesp(d, seed = 8)
  a <- test_accuracy(fixed)
  b <- sensitivity_analysis(mar_fit, odds_ratios = exp(1))
  expect_equal(a$mean, b$mean, tolerance = 0.01)
  # and its intervals contain the true values, while those of MAR do not
  covers <- function(x, q, value) {
    i <- x$quantity == q
    x$lower[i] <= value && value <= x$upper[i]
  }
  a99 <- test_accuracy(fixed, level = 0.99)
  m99 <- test_accuracy(mar_fit, level = 0.99)
  expect_true(covers(a99, "Se[T1]", 0.85) && covers(a99, "prevalence", 0.2))
  expect_false(covers(m99, "prevalence", 0.2))

  # one beta per verified profile
  by_profile <- baysesp(d, saturated(verification = mnar(by_profile = TRUE)), seed = 9)
  expect_identical(sum(startsWith(posterior::variables(by_profile@draws), "beta[")), 4L)
  expect_no_error(sensitivity_analysis(by_profile))
})
