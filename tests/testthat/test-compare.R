# update(), summary() and compare_models().

warnings_of <- function(expr) {
  w <- character()
  withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, conditionMessage(c))
    invokeRestart("muffleWarning")
  })
  w
}

compare_data <- function() {
  profiles <- data.frame(T1 = c(0, 1, 0, 1), T2 = c(0, 0, 1, 1))
  agg <- rbind(
    cbind(profiles, V = 1, D = 1, n = c(1, 8, 5, 30)),
    cbind(profiles, V = 1, D = 0, n = c(49, 22, 25, 15)),
    cbind(profiles, V = 0, D = NA, n = c(750, 70, 20, 5))
  )
  baysesp_data(agg, c("T1", "T2"), "V", "D", counts = "n")
}

test_that("log_marginal() integrates priors on the logit scale", {
  b <- prior_beta(0.3, 5)
  expect_equal(log_marginal(b, 3, 4),
               log(integrate(function(p) p^3 * (1 - p)^4 * dbeta(p, 1.5, 3.5), 0, 1)$value))
  # Riemann sum on a fine grid of the logit scale
  x <- seq(-60, 60, by = 1e-4)
  n <- prior_normal(0.5, 1.2)
  direct <- sum(plogis(x)^30 * plogis(-x)^2 * dnorm(x, 0.5, 1.2)) * 1e-4
  expect_equal(log_marginal(n, 30, 2), log(direct), tolerance = 1e-6)
  s <- prior_student(3, 0, 2)
  direct <- sum(exp(plogis(x, log.p = TRUE) + 200 * plogis(-x, log.p = TRUE)) *
                  dt(x / 2, 3) / 2) * 1e-4
  expect_equal(log_marginal(s, 1, 200), log(direct), tolerance = 1e-6)
  expect_identical(log_marginal(n, 0, 0), 0)
})

test_that("model_cells() leaves out the verified cells of never-verified profiles", {
  cells <- model_cells(example_fit(verified = c(FALSE, rep(TRUE, 7)))@data)
  expect_identical(nrow(cells), 22L)
  expect_identical(cells$cell[1], "unverified")
  expect_identical(cells$label[1], "(0,0,0)")
})

test_that("cell probabilities of the saturated model sum to 1", {
  fit <- example_fit(verified = c(FALSE, rep(TRUE, 7)))
  p <- cell_probs(fit, model_cells(fit@data), prior_beta(0.5, 2))
  expect_equal(rowSums(p), rep(1, nrow(p)))
})

test_that("weighted_total() treats a cell as identical subjects", {
  e <- c(-1, -2, -0.5)
  n <- c(2L, 1L, 3L)
  subjects <- rep(e, n)
  expect_equal(weighted_total(e, n),
               c(sum(subjects), sqrt(length(subjects)) * sd(subjects)))
})

test_that("summary() of a fitted model collects model, priors and diagnostics", {
  fit <- example_fit(verified = c(FALSE, rep(TRUE, 7)))
  fit@model <- saturated(unverified_prevalence_prior = prior_beta(0.1, 10),
                         verification = mnar())
  s <- summary(fit)
  expect_s3_class(s, "fit_summary")
  expect_identical(s$model, "saturated model, transparent parametrization")
  expect_identical(s$verification, "MNAR")
  expect_identical(s$priors$argument,
                   c("frequency_prior", "verified_prior", "profile_prevalence_prior",
                     "unverified_prevalence_prior (0,0,0)", "mnar(odds_ratio_prior)"))
  expect_identical(s$never_verified, "(0,0,0)")
  expect_named(s$diagnostics, c("chains", "draws", "backend", "max_rhat",
                                "min_ess_bulk", "min_ess_tail", "divergent"))
  expect_identical(s$accuracy, rbind(test_accuracy(fit), odds_ratio_summary(fit)))
  expect_output(print(s), "Never verified: \\(0,0,0\\)")
  expect_output(print(s), "unverified_prevalence_prior \\(0,0,0\\) +Beta\\(mean = 0.1, precision = 10\\)")
  expect_error(summary(fit, level = 2), "`level`")
})

test_that("compare_models() checks its arguments", {
  fit <- example_fit()
  expect_error(compare_models(), "one or more models")
  expect_error(compare_models(fit, 1), "one or more models")
  expect_error(compare_models(fit, refit = NA), "`refit`")
  other <- example_fit(verified = c(FALSE, rep(TRUE, 7)))
  expect_error(compare_models(fit, other), "same data")
  expect_error(compare_models(fit, verified_prior = prior_dirichlet(1)), "verified_prior")
})

test_that("posterior_sample() draws from the posterior on the logit scale", {
  set.seed(1)
  x <- posterior_sample(prior_normal(0, 2.5), 30, 10, 1e5)
  grid <- seq(-10, 10, by = 1e-3)
  w <- plogis(grid)^30 * plogis(-grid)^10 * dnorm(grid, 0, 2.5)
  expect_equal(mean(x), sum(plogis(grid) * w) / sum(w), tolerance = 1e-3)
  expect_equal(mean(posterior_sample(prior_beta(0.5, 2), 30, 10, 1e5)), 31 / 42,
               tolerance = 1e-3)
})

test_that("saturated models that differ only in what the data cannot inform are flagged", {
  fit <- example_fit()
  mnar_fit <- fit
  mnar_fit@model <- saturated(verification = mnar())
  # constant draws of a tiny table: WAIC also differs from LOO
  expect_match(warnings_of(compare_models(a = fit, b = mnar_fit)),
               "Models a, b give the same probability", all = FALSE)
  # their own verified_prior is replaced by the common one of compare_models()
  vprior <- fit
  vprior@model <- saturated(verified_prior = prior_beta(0.5, 50))
  set.seed(1)
  w <- character()
  cmp <- withCallingHandlers(compare_models(a = fit, c = vprior), warning = function(c) {
    w <<- c(w, conditionMessage(c))
    invokeRestart("muffleWarning")
  })
  expect_match(w, "Models a, c give the same probability", all = FALSE)
  expect_identical(cmp$table$elpd_loo[1], cmp$table$elpd_loo[2])
  other <- fit
  other@model <- saturated(profile_prevalence_prior = prior_beta(0.5, 10))
  expect_false(any(grepl("same probability", warnings_of(compare_models(fit, other)))))
})

test_that("update() fits a new model with the same settings", {
  skip_on_cran()
  d <- compare_data()
  fit <- baysesp(d, chains = 2, iter_warmup = 300, iter_sampling = 300, seed = 3)
  expect_identical(fit@settings,
                   list(chains = 2, iter_warmup = 300, iter_sampling = 300, seed = 3))
  same <- update(fit)
  expect_identical(same@draws, fit@draws)
  # a short run: only the settings are checked, not the convergence
  ci <- suppressWarnings(update(fit, model = loglinear(), iter_sampling = 200))
  expect_true(S7::S7_inherits(ci@model, loglinear))
  expect_identical(ci@settings$iter_sampling, 200)
  expect_identical(posterior::nchains(ci@draws), 2L)
  expect_identical(posterior::niterations(ci@draws), 200L)
})

test_that("exact LOO of the saturated model agrees with PSIS-LOO", {
  skip_on_cran()
  d <- compare_data()
  fit <- baysesp(d, seed = 4)
  cells <- model_cells(d)
  cells <- cells[cells$count > 0L, ]
  uniform <- prior_beta(0.5, 2)
  exact <- exact_loo(fit, cells, uniform)
  psis <- psis_loo(log(cell_probs(fit, cells, uniform)), 4L)
  expect_true(all(psis$k < 0.7))
  expect_equal(psis$elpd, exact, tolerance = 0.01)
  expect_lt(abs(sum(cells$count * (psis$elpd - exact))), 0.5)

  # a prior on the logit scale: numerical integration
  fit_n <- update(fit, model = saturated(profile_prevalence_prior = prior_normal(0, 1.5)))
  expect_equal(psis_loo(log(cell_probs(fit_n, cells, uniform)), 4L)$elpd,
               exact_loo(fit_n, cells, uniform), tolerance = 0.01)

  # a common verified_prior on the logit scale: grid draws and integration
  logit <- prior_normal(0, 2.5)
  set.seed(4)
  expect_equal(psis_loo(log(cell_probs(fit, cells, logit)), 4L)$elpd,
               exact_loo(fit, cells, logit), tolerance = 0.01)
})

test_that("compare_models() compares structured and saturated models", {
  skip_on_cran()
  set.seed(5)
  d <- compare_data()
  fit <- baysesp(d, seed = 5)
  fit_ci <- update(fit, model = loglinear())
  fit_mnar <- update(fit, model = loglinear(verification = mnar()))
  expect_true(all(paste0("alpha[", 1:4, "]") %in% posterior::variables(fit_mnar@draws)))
  cmp <- compare_models(saturated = fit, fit_ci, fit_mnar)
  expect_s3_class(cmp, "model_comparison")
  expect_setequal(cmp$table$model, c("saturated", "fit_ci", "fit_mnar"))
  expect_identical(cmp$table$elpd_diff[1], 0)
  expect_identical(cmp$table$loo_method[cmp$table$model == "saturated"], "exact")
  expect_true(all(cmp$table$ppc_pvalue >= 0 & cmp$table$ppc_pvalue <= 1))
  expect_identical(nrow(cmp$checks), 3L * 12L)
  expect_true(all(cmp$checks$lower <= cmp$checks$upper))
  expect_output(print(cmp), "Expected log predictive density")
  expect_s3_class(plot(cmp), "ggplot")

  # cell probabilities of the structured models are a distribution
  cells <- model_cells(d)
  uniform <- prior_beta(0.5, 2)
  expect_equal(rowSums(cell_probs(fit_mnar, cells, uniform)), rep(1, 4000))
  expect_equal(rowSums(cell_probs(fit_ci, cells, prior_student(3, 0, 2))), rep(1, 4000))

  # exact refit for one cell agrees with PSIS for a well-behaved cell
  cells <- cells[cells$count > 0L, ]
  psis <- psis_loo(log(cell_probs(fit_ci, cells, uniform)), 4L)
  i <- which(cells$label == "(1,1)" & cells$cell == "diseased")
  expect_equal(refit_loo(fit_ci, cells[i, ], uniform), psis$elpd[i], tolerance = 0.02)
})
