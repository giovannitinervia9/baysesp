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

test_that("dependence = independence() is the default and keeps the old draws", {
  args <- list(n = 200, prevalence = 0.3, sensitivity = c(0.9, 0.7, 0.6),
               specificity = c(0.8, 0.9, 0.95), verification_prob = 0.5)
  set.seed(5)
  d_default <- do.call(simulate_study, args)
  set.seed(5)
  d_explicit <- do.call(simulate_study, c(args, list(dependence = independence())))
  expect_identical(d_default, d_explicit)
  # same draws as the code without the dependence argument
  set.seed(5)
  disease <- stats::rbinom(200, 1L, 0.3)
  p <- outer(disease, args$sensitivity) + outer(1L - disease, 1 - args$specificity)
  results <- matrix(stats::rbinom(600, 1L, p), 200, 3)
  expect_identical(unname(as.matrix(d_default[, 1:3])), results)
})

# theoretical profile probabilities under latent_dependence(), by integration
latent_theta <- function(p, b) {
  K <- length(p)
  a <- stats::qnorm(p) * sqrt(1 + b^2)
  profiles <- as.matrix(expand.grid(rep(list(0:1), K)))
  apply(profiles, 1, function(t) {
    stats::integrate(function(z) {
      vapply(z, function(zz) {
        q <- stats::pnorm(a + b * zz)
        prod(q^t * (1 - q)^(1 - t))
      }, numeric(1)) * stats::dnorm(z)
    }, -Inf, Inf)$value
  })
}

profile_freq <- function(d, tests) {
  j <- 1L + drop(as.matrix(d[, tests]) %*% 2L^(seq_along(tests) - 1L))
  tabulate(j, 2L^length(tests)) / nrow(d)
}

test_that("latent_dependence() gives the theoretical profile probabilities", {
  set.seed(6)
  se <- c(0.9, 0.8, 0.7)
  sp <- c(0.85, 0.9, 0.95)
  b1 <- c(1, 1, -0.5)
  b0 <- c(0.5, 0.5, 0)
  d <- simulate_study(n = 2e5, prevalence = 0.5, sensitivity = se, specificity = sp,
                      verification_prob = 1,
                      dependence = latent_dependence(diseased_loadings = b1,
                                                     healthy_loadings = b0))
  tests <- c("T1", "T2", "T3")
  sick <- d$disease == 1
  expect_equal(profile_freq(d[sick, ], tests), latent_theta(se, b1), tolerance = 0.01)
  # among the healthy P(T = 1) = 1 - Sp
  expect_equal(profile_freq(d[!sick, ], tests), latent_theta(1 - sp, b0), tolerance = 0.01)
  # the dependence is really there: P(T1 = T2 = 1 | D = 1) above Se1 * Se2
  expect_gt(mean(d$T1[sick] * d$T2[sick]), 0.9 * 0.8 + 0.02)
  # marginals unchanged
  expect_equal(unname(colMeans(d[sick, tests])), se, tolerance = 0.01)
  expect_equal(unname(colMeans(1 - d[!sick, tests])), sp, tolerance = 0.01)
})

# P(T_j = 1, T_k = 1) under probit_dependence(): P(Y_j < c_j, Y_k < c_k)
both_positive <- function(pj, pk, r) {
  cj <- stats::qnorm(pj)
  ck <- stats::qnorm(pk)
  stats::integrate(function(y) {
    stats::dnorm(y) * stats::pnorm((ck - r * y) / sqrt(1 - r^2))
  }, -Inf, cj)$value
}

test_that("probit_dependence() gives the theoretical pairwise probabilities", {
  set.seed(7)
  se <- c(0.9, 0.8, 0.7)
  sp <- c(0.85, 0.9, 0.95)
  R1 <- matrix(c(1, 0.5, -0.3, 0.5, 1, -0.3, -0.3, -0.3, 1), 3, 3)
  d <- simulate_study(n = 2e5, prevalence = 0.5, sensitivity = se, specificity = sp,
                      verification_prob = 1,
                      dependence = probit_dependence(diseased_correlation = R1,
                                                     healthy_correlation = -0.3))
  sick <- d$disease == 1
  x1 <- as.matrix(d[sick, 1:3])
  x0 <- as.matrix(d[!sick, 1:3])
  for (pair in list(c(1, 2), c(1, 3), c(2, 3))) {
    j <- pair[1]
    k <- pair[2]
    expect_equal(mean(x1[, j] * x1[, k]), both_positive(se[j], se[k], R1[j, k]),
                 tolerance = 0.01)
    expect_equal(mean(x0[, j] * x0[, k]), both_positive(1 - sp[j], 1 - sp[k], -0.3),
                 tolerance = 0.01)
  }
  expect_equal(unname(colMeans(x1)), se, tolerance = 0.01)
  expect_equal(unname(colMeans(1 - x0)), sp, tolerance = 0.01)
})

test_that("dependence mechanisms reject invalid values", {
  expect_error(latent_dependence(diseased_loadings = NA_real_), "finite")
  expect_error(latent_dependence(healthy_loadings = Inf), "finite")
  expect_error(probit_dependence(diseased_correlation = 1.5), "[-1, 1]", fixed = TRUE)
  expect_error(probit_dependence(diseased_correlation = c(0.1, 0.2)), "single number")
  expect_error(probit_dependence(healthy_correlation = matrix(c(1, 0.2, 0.3, 1), 2)),
               "symmetric")
  args <- list(n = 10, prevalence = 0.3, sensitivity = c(0.9, 0.7, 0.6),
               specificity = c(0.8, 0.9, 0.95), verification_prob = 1)
  expect_error(do.call(simulate_study, c(args, list(dependence = "none"))),
               "`dependence`")
  expect_error(do.call(simulate_study, c(args, list(
    dependence = latent_dependence(diseased_loadings = c(1, 1))))), "diseased_loadings")
  expect_error(do.call(simulate_study, c(args, list(
    dependence = probit_dependence(diseased_correlation = -0.6)))), "positive definite")
  expect_error(do.call(simulate_study, c(args, list(
    dependence = probit_dependence(healthy_correlation = diag(2))))), "3 x 3")
})
