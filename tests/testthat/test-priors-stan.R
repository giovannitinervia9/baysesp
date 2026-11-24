# Checks that prior_lpdf() in inst/stan/include/priors.stan defines the intended
# densities: a parameter sampled from the prior alone must have the moments of
# the prior. Needs cmdstanr and CmdStan, so it does not run on CRAN.

sample_prior <- function(prior) {
  code <- prior_code(prior)
  model_file <- tempfile(fileext = ".stan")
  writeLines(c(
    "functions {",
    "#include priors.stan",
    "}",
    "data { int family; vector[3] params; }",
    "parameters { real x; }",
    "model { target += prior_lpdf(x | family, params); }"
  ), model_file)
  model <- cmdstanr::cmdstan_model(
    model_file,
    include_paths = system.file("stan", "include", package = "baysesp"),
    quiet = TRUE
  )
  fit <- model$sample(
    data = list(family = code$family, params = code$params),
    chains = 4, iter_warmup = 500, iter_sampling = 5000,
    seed = 2026, refresh = 0, show_messages = FALSE
  )
  as.numeric(fit$draws("x", format = "draws_matrix"))
}

test_that("prior_lpdf reproduces the normal, Student t and beta densities", {
  skip_on_cran()
  skip_if_not(baysesp:::cmdstan_available(), "cmdstanr and CmdStan are not available")

  x <- sample_prior(prior_normal(1, 2))
  expect_equal(mean(x), 1, tolerance = 0.05)
  expect_equal(sd(x), 2, tolerance = 0.05)

  x <- sample_prior(prior_student(df = 5, location = -1, scale = 0.5))
  expect_equal(median(x), -1, tolerance = 0.05)
  expect_equal(sd(x), 0.5 * sqrt(5 / 3), tolerance = 0.05)

  # beta: the density is on p = inv_logit(x), with the Jacobian
  p <- stats::plogis(sample_prior(prior_beta(mean = 0.2, precision = 10)))
  expect_equal(mean(p), 0.2, tolerance = 0.03)
  expect_equal(var(p), 0.2 * 0.8 / 11, tolerance = 0.05)
})
