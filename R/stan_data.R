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

# In the saturated model beta only links the proportion of diseased among the
# unverified to that among the verified (specification, Proposition 7.2).
S7::method(stan_data, list(saturated, mnar)) <- function(model, verification, data) {
  s <- saturated_data(model, data)
  beta <- prior_code(verification@odds_ratio_prior)
  s$n_beta <- if (verification@by_profile) s$J - s$J0 else 1L
  s$beta_family <- beta$family
  s$beta_params <- beta$params
  s
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
  q <- lapply(q_priors(q, J0, "unverified_prevalence_prior"), prior_code)
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
    n_beta = 0L,
    beta_family = 1L,
    beta_params = c(0, 1, 0)
  )
}

# One prior, or a list with one prior per never-verified profile, as a list of
# length J0.
q_priors <- function(q, J0, argument) {
  if (length(q) > 1L && length(q) != J0) {
    stop("`", argument, "` must be one prior, or a list with one prior ",
         "per never-verified profile (", J0, ").", call. = FALSE)
  }
  rep_len(q, J0)
}

S7::method(stan_program, loglinear) <- function(model) "loglinear"

S7::method(stan_data, list(loglinear, mar)) <- function(model, verification, data) {
  loglinear_data(model, data)
}

S7::method(stan_data, list(loglinear, mnar)) <- function(model, verification, data) {
  structured_mnar(loglinear_data(model, data), verification, data)
}

# In the structured models beta enters the verification probabilities, with
# one intercept per verified profile (specification, section 6.1). The data
# lists of loglinear.stan and random_effects.stan share these elements.
structured_mnar <- function(s, verification, data) {
  alpha <- prior_code(verification@intercepts_prior)
  beta <- prior_code(verification@odds_ratio_prior)
  s$mnar <- 1L
  s$n_beta <- if (verification@by_profile) sum(!data@never_verified) else 1L
  s$alpha_family <- alpha$family
  s$alpha_params <- alpha$params
  s$beta_family <- beta$family
  s$beta_params <- beta$params
  s
}

structured_mar <- list(mnar = 0L, n_beta = 0L, beta_value = 0,
                       alpha_family = 1L, alpha_params = c(0, 1, 0),
                       beta_family = 1L, beta_params = c(0, 1, 0))

# Data list of loglinear.stan, MAR.
loglinear_data <- function(model, data) {
  K <- length(data@tests)
  if (model@interaction_order > K) {
    stop("`interaction_order` (", model@interaction_order, ") cannot exceed the ",
         "number of tests (", K, ").", call. = FALSE)
  }
  X <- design_matrix(data@profiles, model@interaction_order)
  prevalence <- prior_code(model@prevalence_prior)
  main <- prior_code(model@main_effects_prior)
  interaction <- prior_code(model@interactions_prior)
  c(list(
    J = nrow(X), P = ncol(X), X = unname(X),
    order = as.array(attr(X, "order")),
    n1 = data@n_diseased, n0 = data@n_healthy, u = data@n_unverified,
    prevalence_family = prevalence$family, prevalence_params = prevalence$params,
    main_family = main$family, main_params = main$params,
    interaction_family = interaction$family, interaction_params = interaction$params
  ), structured_mar)
}

# Design matrix X_m (specification, section 5.2): one column per subset S of at
# most m tests, by size, with x_S(t) = prod over k in S of t_k. Columns are
# named "A", "A:B", ...; attribute "order" is the size of each subset.
design_matrix <- function(profiles, m) {
  K <- ncol(profiles)
  subsets <- unlist(lapply(seq_len(m), function(s) utils::combn(K, s, simplify = FALSE)),
                    recursive = FALSE)
  X <- vapply(subsets, function(S) as.numeric(apply(profiles[, S, drop = FALSE], 1, prod)),
              numeric(nrow(profiles)))
  X <- matrix(X, nrow = nrow(profiles))
  colnames(X) <- vapply(subsets, function(S) paste(colnames(profiles)[S], collapse = ":"),
                        character(1))
  attr(X, "order") <- lengths(subsets)
  X
}

S7::method(stan_program, random_effects) <- function(model) "random_effects"

S7::method(stan_data, list(random_effects, mar)) <- function(model, verification, data) {
  random_effects_data(model, data)
}

S7::method(stan_data, list(random_effects, mnar)) <- function(model, verification, data) {
  structured_mnar(random_effects_data(model, data), verification, data)
}

# Data list of random_effects.stan, MAR.
random_effects_data <- function(model, data) {
  K <- length(data@tests)
  if (K < 2L) {
    stop("The random effects model needs at least two tests.", call. = FALSE)
  }
  gh <- gauss_hermite(32L)
  prevalence <- prior_code(model@prevalence_prior)
  intercepts <- prior_code(model@intercepts_prior)
  loadings <- prior_code(model@loadings_prior)
  c(list(
    J = nrow(data@profiles), K = K, profiles = unname(data@profiles * 1),
    positive_loadings = as.integer(model@positive_loadings),
    Q = length(gh$nodes), nodes = gh$nodes, log_weights = log(gh$weights),
    n1 = data@n_diseased, n0 = data@n_healthy, u = data@n_unverified,
    prevalence_family = prevalence$family, prevalence_params = prevalence$params,
    intercepts_family = intercepts$family, intercepts_params = intercepts$params,
    loadings_family = loadings$family, loadings_params = loadings$params
  ), structured_mar)
}

# Gauss-Hermite quadrature for the standard normal distribution (Golub-Welsch):
# the nodes are the eigenvalues of the Jacobi matrix of the probabilists'
# Hermite polynomials, the weights the squared first components of the
# eigenvectors. sum(weights * f(nodes)) approximates E f(Z), Z ~ N(0, 1).
gauss_hermite <- function(n) {
  jacobi <- matrix(0, n, n)
  off <- sqrt(seq_len(n - 1L))
  jacobi[cbind(1:(n - 1L), 2:n)] <- off
  jacobi[cbind(2:n, 1:(n - 1L))] <- off
  e <- eigen(jacobi, symmetric = TRUE)
  list(nodes = rev(e$values), weights = rev(e$vectors[1, ]^2))
}
