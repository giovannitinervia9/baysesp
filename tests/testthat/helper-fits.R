# A baysesp_fit with constant draws, for testing the functions used after the
# fit without running Stan. Profile distributions of Example 4.5 of the
# specification (K = 3 tests A, B, C), prevalence 0.1. `verified` says which
# profiles have verified subjects. The transparent parameters are consistent
# with MAR: phi1 = psi, with half of each verified profile verified.
example_theta1 <- c(0.05, 0.08, 0.04, 0.10, 0.03, 0.12, 0.08, 0.50)
example_theta0 <- c(0.80, 0.06, 0.05, 0.02, 0.04, 0.01, 0.01, 0.01)

example_fit <- function(verified = rep(TRUE, 8), prevalence = 0.1) {
  omega <- prevalence * example_theta1 + (1 - prevalence) * example_theta0
  psi <- prevalence * example_theta1 / omega
  named <- function(x, name) {
    if (length(x) == 0L) return(list())
    stats::setNames(as.list(x), paste0(name, "[", seq_along(x), "]"))
  }
  vars <- c(list(prevalence = prevalence),
            named(example_theta0, "theta_healthy"), named(example_theta1, "theta_diseased"),
            named(omega, "omega"), named(ifelse(verified, 0.5, 0), "v"),
            named(ifelse(verified, psi, 0), "phi1"), named(psi[!verified], "q"))
  draws <- do.call(posterior::draws_array, lapply(vars, rep, 4))
  x <- data.frame(expand.grid(A = 0:1, B = 0:1, C = 0:1), V = as.integer(verified), D = 1)
  d <- baysesp_data(x, c("A", "B", "C"), "V", "D")
  baysesp_fit(model = saturated(), data = d, draws = draws, backend = "rstan",
              stan_fit = NULL)
}
