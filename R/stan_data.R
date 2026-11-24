#' @include models.R data.R
NULL

# From model and data to the Stan programs ---------------------------------------
#
# stan_data() dispatches on both the model and the verification mechanism,
# because the same assumption is translated differently in the two roads
# (specification, section 7.5). stan_program() gives the name of the program.

stan_data <- S7::new_generic("stan_data", c("model", "verification"),
                             function(model, verification, data) S7::S7_dispatch())

stan_program <- S7::new_generic("stan_program", "model")

S7::method(stan_program, saturated) <- function(model) "saturated"

S7::method(stan_data, list(saturated, mar)) <- function(model, verification, data) {
  saturated_data(model, data)
}

# Data list of saturated.stan, MAR; mnar() changes only the last three elements.
saturated_data <- function(model, data) {
  J <- length(data@n_diseased)
  never <- data@never_verified
  J0 <- sum(never)

  concentration <- model@frequency_prior@concentration
  if (!length(concentration) %in% c(1L, J)) {
    stop("`frequency_prior` must have one concentration, or one per profile (", J, ").",
         call. = FALSE)
  }

  q <- as_prior_list(model@unverified_prevalence_prior, "unverified_prevalence_prior")
  if (J0 > 0L && length(q) == 0L) {
    stop("No subject is verified in the profiles ",
         paste(profile_labels(data)[never], collapse = " "), ". ",
         "The proportion of diseased in these profiles cannot be estimated from the ",
         "data: give it a prior with `saturated(unverified_prevalence_prior = )`, ",
         "for example with elicit_beta().", call. = FALSE)
  }
  if (length(q) > 1L && length(q) != J0) {
    stop("`unverified_prevalence_prior` must be one prior, or a list with one prior ",
         "per never-verified profile (", J0, ").", call. = FALSE)
  }
  q <- lapply(rep_len(q, J0), prior_code)
  verified <- prior_code(model@verified_prior)
  prevalence <- prior_code(model@profile_prevalence_prior)

  list(
    J = J,
    n1 = data@n_diseased,
    n0 = data@n_healthy,
    u = data@n_unverified,
    frequency_concentration = rep_len(concentration, J),
    verified_family = verified$family,
    verified_params = verified$params,
    prevalence_family = prevalence$family,
    prevalence_params = prevalence$params,
    J0 = J0,
    # as.array(): rstan needs arrays of length 1 to be arrays
    q_family = as.array(vapply(q, `[[`, integer(1), "family")),
    q_params = matrix(as.double(unlist(lapply(q, `[[`, "params"))), J0, 3, byrow = TRUE),
    mnar = 0L,
    beta_family = 1L,
    beta_params = c(0, 1, 0)
  )
}
