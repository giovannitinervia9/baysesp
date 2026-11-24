test_that("the Stan programs are shipped with the package", {
  expect_true("saturated" %in% baysesp:::stan_programs())
  expect_true(file.exists(baysesp:::stan_file("saturated")))
  expect_error(baysesp:::stan_file("does_not_exist"), "not found")
})

test_that("rstan is the default backend", {
  expect_identical(baysesp:::resolve_backend(), "rstan")
  expect_true(suppressMessages(baysesp:::resolve_backend("cmdstanr")) %in% c("cmdstanr", "rstan"))
})

test_that("the rstan backend has the precompiled programs", {
  expect_true("saturated" %in% names(baysesp:::stanmodels))
})

test_that("the Stan hash depends only on the code", {
  h <- baysesp:::stan_hash("saturated")
  expect_type(h, "character")
  expect_identical(nchar(h), 12L)
  expect_identical(h, baysesp:::stan_hash("saturated"))
})
