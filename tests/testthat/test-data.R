# aggregated table of Example 2.1 of the specification (K = 2, N = 1000)
example_table <- function() {
  profiles <- data.frame(T1 = c(0, 1, 0, 1), T2 = c(0, 0, 1, 1))
  rbind(
    cbind(profiles, V = 1, D = 1, n = c(1, 8, 5, 30)),
    cbind(profiles, V = 1, D = 0, n = c(49, 22, 25, 15)),
    cbind(profiles, V = 0, D = NA, n = c(750, 70, 20, 5))
  )
}

test_that("baysesp_data counts the cells of aggregated data", {
  d <- baysesp_data(example_table(), tests = c("T1", "T2"), verified = "V",
                    disease = "D", counts = "n")
  expect_true(S7::S7_inherits(d, baysesp_data))
  expect_identical(d@tests, c("T1", "T2"))
  expect_identical(d@profiles,
                   matrix(c(0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L), 4, 2,
                          dimnames = list(NULL, c("T1", "T2"))))
  expect_identical(d@n_diseased, c(1L, 8L, 5L, 30L))
  expect_identical(d@n_healthy, c(49L, 22L, 25L, 15L))
  expect_identical(d@n_unverified, c(750L, 70L, 20L, 5L))
  expect_identical(d@never_verified, rep(FALSE, 4))
})

test_that("individual and aggregated data give the same object", {
  agg <- example_table()
  ind <- agg[rep(seq_len(nrow(agg)), agg$n), c("T1", "T2", "V", "D")]
  # logical coding and a non-missing disease for unverified subjects are accepted
  ind$T2 <- ind$T2 == 1
  ind$D[ind$V == 0] <- 1
  expect_identical(
    baysesp_data(ind, tests = c("T1", "T2"), verified = "V", disease = "D"),
    baysesp_data(agg, tests = c("T1", "T2"), verified = "V", disease = "D",
                 counts = "n")
  )
})

test_that("profiles are in binary order with the first test fastest", {
  x <- data.frame(A = 1, B = 0, C = 1, V = 1, D = 1)
  d <- baysesp_data(x, tests = c("A", "B", "C"), verified = "V", disease = "D")
  expect_identical(nrow(d@profiles), 8L)
  # profile (1, 0, 1) is in row 1 + 1 + 4 = 6
  expect_identical(which(d@n_diseased == 1L), 6L)
  expect_identical(unname(d@profiles[6, ]), c(1L, 0L, 1L))
  expect_identical(sum(d@never_verified), 7L)
})

test_that("a single test works", {
  x <- data.frame(T = c(0, 1, 1, 0), V = c(1, 1, 0, 0), D = c(0, 1, NA, NA))
  d <- baysesp_data(x, tests = "T", verified = "V", disease = "D")
  expect_identical(d@n_diseased, c(0L, 1L))
  expect_identical(d@n_healthy, c(1L, 0L))
  expect_identical(d@n_unverified, c(1L, 1L))
  expect_output(print(d), "1 test, 2 profiles, 4 subjects, 2 verified")
})

test_that("never-verified profiles are detected and printed", {
  x <- data.frame(T1 = c(1, 0, 1), T2 = c(1, 0, 0), V = c(1, 0, 1), D = c(1, NA, 0))
  d <- baysesp_data(x, tests = c("T1", "T2"), verified = "V", disease = "D")
  expect_identical(d@never_verified, c(TRUE, FALSE, TRUE, FALSE))
  expect_output(print(d), "never-verified profiles: \\(0,0\\) \\(0,1\\)")
  expect_identical(names(as.data.frame(d)),
                   c("T1", "T2", "diseased", "healthy", "unverified"))
})

test_that("baysesp_data rejects invalid data", {
  x <- data.frame(T1 = c(0, 1), V = c(1, 0), D = c(1, NA))
  expect_error(baysesp_data(as.list(x), "T1", "V", "D"), "`data`")
  expect_error(baysesp_data(x, "T9", "V", "D"), "`T9`")
  expect_error(baysesp_data(x, c("T1", "T1"), "V", "D"), "`tests`")
  expect_error(baysesp_data(transform(x, T1 = c(0, 2)), "T1", "V", "D"), "0/1")
  expect_error(baysesp_data(transform(x, T1 = c("a", "b")), "T1", "V", "D"), "0/1")
  expect_error(baysesp_data(transform(x, T1 = c(0, NA)), "T1", "V", "D"), "missing")
  expect_error(baysesp_data(transform(x, V = c(1, NA)), "T1", "V", "D"), "missing")
  expect_error(baysesp_data(transform(x, D = c(NA, NA)), "T1", "V", "D"),
               "missing for verified")
  expect_error(baysesp_data(transform(x, V = c(0, 0)), "T1", "V", "D"),
               "at least one subject must be verified")
  expect_error(baysesp_data(transform(x, n = c(1.5, 2)), "T1", "V", "D", counts = "n"),
               "whole numbers")
  expect_error(baysesp_data(transform(x, n = c(-1, 2)), "T1", "V", "D", counts = "n"),
               "whole numbers")
})
