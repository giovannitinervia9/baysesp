#' @include models.R data.R generics.R
NULL

#' Bayesian estimation of the accuracy of diagnostic tests
#'
#' Fits a Bayesian model for the sensitivity and specificity of one or more
#' binary tests when not all subjects have been verified with the reference
#' standard. The objects of the model classes passed to `model` contain the
#' specification of the priors and of the verification mechanism.
#'
#' @details
#' The posterior distribution of the parameters of `model`, given `data`, is
#' approximated with the No-U-Turn sampler of Stan, a form of Hamiltonian Monte
#' Carlo that adapts its step size and mass matrix during the warmup iterations.
#' The likelihood and the priors are coded in a Stan program specific to the
#' model class and to the verification mechanism. Each chain runs
#' `iter_warmup` warmup iterations, which are discarded, followed by
#' `iter_sampling` sampling iterations, which are kept as draws. The draws
#' therefore form `chains * iter_sampling` joint samples of the parameters of
#' the model, and every quantity reported by the package is computed from them.
#'
#' @param data A `baysesp_data` object, created by [baysesp_data()], with the
#'   results of the tests, the verification indicator and the disease status of
#'   each subject. Any other object gives an error.
#' @param model An accuracy model: an object created by one of the model
#'   constructors [saturated()], [loglinear()] or [random_effects()]. It holds
#'   the priors on the model parameters and the verification mechanism, which is
#'   `verification = mar()` or `verification = mnar()` given to the constructor. Default
#'   `saturated()`, with MAR verification. Any other object gives an error.
#' @param backend Stan interface, a single value: `"rstan"` (default) or
#'   `"cmdstanr"`. With `"rstan"`, the Stan programs are the ones compiled when
#'   the package is installed, so nothing is compiled at run time. With
#'   `"cmdstanr"`, each program is compiled the first time it is used and the
#'   executable is kept in the user cache directory, so the first run is slower
#'   and later runs reuse the compiled program. If `cmdstanr` or CmdStan is not
#'   available, `"rstan"` is used and a message is printed. The `backend` element
#'   of the result reports the backend that was actually used. Any other value
#'   gives an error.
#' @param chains Positive integer, the number of Markov chains, which run
#'   independently. Default `4`.
#' @param iter_warmup Positive integer, the number of warmup iterations of each
#'   chain. The warmup iterations adapt the sampler and are then discarded.
#'   Default `1000`.
#' @param iter_sampling Positive integer, the number of sampling iterations of
#'   each chain, kept as draws. Default `1000`, which gives `chains * 1000`
#'   draws in total.
#' @param seed `NULL` (default) or a single whole number between `0` and
#'   `2147483647`. If `NULL`, a seed is drawn from the random number generator
#'   of R, so that calling `set.seed()` before `baysesp()` makes the run
#'   reproducible. A number is used directly as the seed of the sampler.
#' @param ... Further arguments passed to [rstan::sampling()] when `backend` is
#'   `"rstan"` (for example `cores`, to run the chains in parallel, or
#'   `control`, for example `list(adapt_delta = 0.95)`), or to the `$sample()`
#'   method of `cmdstanr` when `backend` is `"cmdstanr"` (for example
#'   `parallel_chains`, or `adapt_delta`). Do not repeat `chains`,
#'   `iter_warmup`, `iter_sampling` or `seed`, which `baysesp()` already passes to
#'   the backend.
#'
#' @return An object of class `baysesp_fit` with the following properties:
#'
#'   * `model`: the model object given in `model`.
#'   * `data`: the `baysesp_data` object given in `data`.
#'   * `draws`: a `posterior::draws_array` with the sampling iterations only,
#'     that is `iter_sampling` iterations, `chains` chains and the variables
#'     listed below.
#'   * `backend`: the backend used, `"rstan"` or `"cmdstanr"`.
#'   * `stan_fit`: the object returned by the backend: a `stanfit` object for
#'     `"rstan"`, a `CmdStanMCMC` object for `"cmdstanr"`.
#'
#' The profiles are the \eqn{J = 2^K} response profiles of the \eqn{K} tests,
#' indexed in binary order: profile \eqn{j} is the pattern of results with index
#' \eqn{j = 1 + \sum_k t_k 2^{k-1}}, where \eqn{t_k = 1} for a positive result on
#' test \eqn{k} and \eqn{t_k = 0} otherwise. Profile 1 has all tests negative,
#' and the first test changes fastest. The variables in `draws` are:
#'
#' * `prevalence`: the prevalence of the disease.
#' * `theta_healthy[j]` and `theta_diseased[j]`: the probability of profile
#'   \eqn{j} among the healthy and among the diseased. Within each draw, each of
#'   the two sets sums to 1 over \eqn{j}.
#' * `omega[j]`, `v[j]` and `phi1[j]`, only for the saturated model: the
#'   probability of profile \eqn{j}, the probability of verification given
#'   profile \eqn{j}, and the probability of disease among the verified subjects
#'   with profile \eqn{j}. `v[j]` and `phi1[j]` are 0 for profiles with no
#'   verified subjects.
#' * `q[b]`: the prevalence of disease among the subjects with the \eqn{b}-th
#'   profile that has no verified subjects. Present only if some profile has no
#'   verified subjects.
#' * `beta`: the log odds ratio of verification between diseased and healthy
#'   subjects with the same profile, present only with MNAR verification. Without
#'   `by_profile = TRUE` it is a single value. With `by_profile = TRUE` there is
#'   one value, named `beta[i]`, for each profile with verified subjects, in
#'   binary order.
#' @seealso [test_accuracy()]
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- simulate_study(n = 1000, prevalence = 0.2,
#'                     sensitivity = c(A = 0.9, B = 0.8),
#'                     specificity = c(A = 0.85, B = 0.95),
#'                     verification_prob = c(0.2, 0.9, 0.9, 1))
#' d <- baysesp_data(d, tests = c("A", "B"), verified = "verified", disease = "disease")
#' fit <- baysesp(d, seed = 1)
#' test_accuracy(fit)
#' }
#' @export
baysesp <- function(data, model = saturated(), backend = c("rstan", "cmdstanr"),
                    chains = 4, iter_warmup = 1000, iter_sampling = 1000,
                    seed = NULL, ...) {
  if (!S7::S7_inherits(data, baysesp_data)) {
    stop("`data` must be created with baysesp_data().", call. = FALSE)
  }
  if (!S7::S7_inherits(model, accuracy_model)) {
    stop("`model` must be a model, such as saturated().", call. = FALSE)
  }
  backend <- resolve_backend(backend)
  sdata <- stan_data(model, model@verification, data)
  fit_stan(model, data, sdata, backend, chains, iter_warmup, iter_sampling, seed, ...)
}

# Runs the Stan program of the model on a data list and returns a baysesp_fit.
fit_stan <- function(model, data, sdata, backend, chains, iter_warmup, iter_sampling,
                     seed, ...) {
  program <- stan_program(model)
  variables <- c("prevalence", "theta_healthy", "theta_diseased",
                 if (S7::S7_inherits(model, saturated)) c("omega", "v", "phi1"),
                 if (isTRUE(sdata$J0 > 0L)) "q", if (sdata$n_beta > 0L) "beta")

  if (backend == "rstan") {
    if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
    stan_fit <- rstan::sampling(
      stanmodels[[program]], data = sdata, pars = variables,
      chains = chains, iter = iter_warmup + iter_sampling, warmup = iter_warmup,
      seed = seed, refresh = 0, ...
    )
    draws <- posterior::as_draws_array(rstan::extract(stan_fit, pars = variables,
                                                      permuted = FALSE))
  } else {
    stan_fit <- cmdstan_program(program)$sample(
      data = sdata, chains = chains, iter_warmup = iter_warmup,
      iter_sampling = iter_sampling, seed = seed, refresh = 0,
      show_messages = FALSE, ...
    )
    draws <- stan_fit$draws(variables)
  }

  baysesp_fit(model = model, data = data, draws = draws, backend = backend,
              stan_fit = stan_fit)
}

#' @rdname baysesp
#' @usage NULL
#' @export
baysesp_fit <- S7::new_class(
  "baysesp_fit",
  properties = list(
    model = accuracy_model,
    data = baysesp_data,
    draws = S7::class_any,
    backend = S7::class_character,
    stan_fit = S7::class_any
  )
)

S7::method(print_generic, baysesp_fit) <- function(x, ...) {
  cat("<baysesp_fit> ", sub("^baysesp::", "", class(x@model)[1]), " model, ",
      toupper(sub("^baysesp::", "", class(x@model@verification)[1])),
      " verification, ", posterior::nchains(x@draws), " chains x ",
      posterior::niterations(x@draws), " draws (", x@backend, ")\n", sep = "")
  cat(convergence_line(x), "\n\n", sep = "")
  print(rbind(test_accuracy(x), odds_ratio_summary(x)), row.names = FALSE, digits = 3)
  if (S7::S7_inherits(x@model@verification, mar)) {
    cat("\nOdds ratio of verification fixed at 1 (MAR).\n")
  }
  invisible(x)
}

# One line with the worst R-hat and ESS over all saved variables and the number
# of divergent transitions; constant variables (R-hat NA) are skipped.
convergence_line <- function(fit) {
  d <- posterior::summarise_draws(fit@draws, "rhat", "ess_bulk", "ess_tail")
  worst <- function(v, fun, digits) {
    v <- v[is.finite(v)]
    if (length(v) == 0L) "NA" else formatC(fun(v), digits = digits, format = "f")
  }
  sprintf("Convergence: max R-hat %s, min bulk ESS %s, min tail ESS %s, %s divergent transitions",
          worst(d$rhat, max, 3), worst(d$ess_bulk, min, 0), worst(d$ess_tail, min, 0),
          divergent_transitions(fit))
}

divergent_transitions <- function(fit) {
  s <- fit@stan_fit
  if (inherits(s, "stanfit")) return(rstan::get_num_divergent(s))
  if (inherits(s, "CmdStanMCMC")) {
    return(sum(s$diagnostic_summary("divergences", quiet = TRUE)$num_divergent))
  }
  NA_integer_
}

# Posterior summary of the odds ratio of verification (exp(beta)), estimated
# only with MNAR verification and never identified by the data.
odds_ratio_summary <- function(fit, level = 0.95) {
  if (!any(grepl("^beta(\\[|$)", posterior::variables(fit@draws)))) return(NULL)
  b <- unclass(posterior::as_draws_matrix(posterior::subset_draws(fit@draws, "beta")))
  or <- exp(b)
  colnames(or) <- if (ncol(or) == 1L) "odds_ratio" else sub("^beta", "odds_ratio", colnames(b))
  summarise_draws(or, level, identified = FALSE)
}

# Prevalence and distributions of the profiles, the parameters shared by all
# models, as a vector and two draws x profiles matrices.
core_draws <- function(fit) {
  m <- unclass(posterior::as_draws_matrix(fit@draws))
  J <- nrow(fit@data@profiles)
  list(
    prevalence = unname(m[, "prevalence"]),
    theta_healthy = unname(m[, paste0("theta_healthy[", seq_len(J), "]"), drop = FALSE]),
    theta_diseased = unname(m[, paste0("theta_diseased[", seq_len(J), "]"), drop = FALSE])
  )
}
