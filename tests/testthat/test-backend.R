test_that("the Stan programs are shipped with the package", {
  expect_true("saturated" %in% baysesp:::stan_programs())
  expect_true(file.exists(baysesp:::stan_file("saturated")))
  expect_error(baysesp:::stan_file("non_esiste"), "not found")
})

test_that("the backend is resolved to a valid choice", {
  expect_true(baysesp:::resolve_backend("rstan") == "rstan")
  expect_true(suppressMessages(baysesp:::resolve_backend()) %in% c("cmdstanr", "rstan"))
})

test_that("the rstan backend has the precompiled programs", {
  expect_true("saturated" %in% names(baysesp:::stanmodels))
})
