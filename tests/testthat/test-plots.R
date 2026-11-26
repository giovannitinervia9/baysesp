test_that("plots of priors and fits run", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  for (p in list(plot(prior_beta(0.2, 10)), plot(prior_normal(-1, 0.5)),
                 plot(prior_student(3, 0, 1)), plot(example_fit()))) {
    expect_true(ggplot2::is_ggplot(p))
    expect_no_error(print(p))
  }
  expect_error(plot(prior_dirichlet(1)), "marginals")
})

test_that("prior_posterior() shows that q is not learned from the data", {
  skip_on_cran()
  set.seed(51)
  s <- simulate_study(n = 2000, prevalence = 0.2, sensitivity = c(0.9, 0.7),
                      specificity = c(0.85, 0.95), verification_prob = c(0, 0.5, 0.5, 1))
  d <- baysesp_data(s, c("T1", "T2"), "verified", "disease")
  fit <- baysesp(d, saturated(unverified_prevalence_prior = prior_beta(0.05, 50)), seed = 52)
  pp <- prior_posterior(fit, seed = 53)
  expect_s3_class(pp, "prior_posterior")
  expect_identical(names(pp), c("quantity", "distribution", "mean", "median", "lower", "upper"))
  expect_identical(pp$distribution[1:2], c("prior", "posterior"))
  # q: posterior equal to the prior; PPV[T1]: learned from the data
  q <- pp[pp$quantity == "q[1]", ]
  expect_equal(q$mean[1], q$mean[2], tolerance = 0.05)
  ppv <- pp[pp$quantity == "PPV[T1]", ]
  expect_lt(ppv$upper[2] - ppv$lower[2], (ppv$upper[1] - ppv$lower[1]) / 3)
  expect_output(print(pp), "posterior")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_true(ggplot2::is_ggplot(plot(pp)))
  expect_no_error(print(plot(pp)))
})
