#' @include models.R data.R generics.R
NULL

#' Bayesian estimation of the accuracy of diagnostic tests
#'
#' Fits a Bayesian model for the sensitivity and specificity of one or more
#' binary tests when not all subjects have been verified with the reference
#' standard (Baker, 1995; Buzoianu and Kadane, 2008; Martinez et al., 2006). The objects of the model classes passed to `model` contain the
#' specification of the priors and of the verification mechanism.
#'
#' @details
#' The posterior distribution of the parameters of `model`, given `data`, is
#' approximated with the No-U-Turn sampler of Stan, a form of Hamiltonian Monte
#' Carlo that adapts its step size and mass matrix during the warmup iterations (Carpenter et al., 2017).
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
#'   `"cmdstanr"`. With `"rstan"` (Stan Development Team, 2025), the Stan programs are the ones compiled when
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
#'   `parallel_chains`, or `adapt_delta`). With [saturated()] the target
#'   acceptance rate `adapt_delta` is 0.95 unless given here, instead of the
#'   Stan default 0.8: see the `profile_prevalence_prior` argument of
#'   [saturated()]. Do not repeat `chains`,
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
#'   * `settings`: a list with `chains`, `iter_warmup`, `iter_sampling`, the
#'     seed actually used and the arguments given in `...`. `update()` reuses
#'     them to fit a different model to the same data.
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
#'   subjects with the same profile, present only with MNAR verification (Kosinski and Barnhart, 2003). Without
#'   `by_profile = TRUE` it is a single value. With `by_profile = TRUE` there is
#'   one value, named `beta[i]`, for each profile with verified subjects, in
#'   binary order.
#' * `alpha[i]`, only for [loglinear()] and [random_effects()] with MNAR
#'   verification: the log odds of verification of the healthy subjects with
#'   the \eqn{i}-th profile that has verified subjects, in binary order.
#' @references
#' Baker, S. G. (1995). Evaluating multiple diagnostic tests with partial verification. *Biometrics*, 51(1), 330–337. \doi{10.2307/2533339}
#'
#' Buzoianu, M. and Kadane, J. B. (2008). Adjusting for verification bias in diagnostic test evaluation: A Bayesian approach. *Statistics in Medicine*, 27(13), 2453–2473. \doi{10.1002/sim.3099}
#'
#' Carpenter, B., Gelman, A., Hoffman, M. D., Lee, D., Goodrich, B., Betancourt, M., Brubaker, M., Guo, J., Li, P. and Riddell, A. (2017). Stan: A probabilistic programming language. *Journal of Statistical Software*, 76(1), 1–32. \doi{10.18637/jss.v076.i01}
#'
#' Kosinski, A. S. and Barnhart, H. X. (2003). Accounting for nonignorable verification bias in assessment of diagnostic tests. *Biometrics*, 59(1), 163–171. \doi{10.1111/1541-0420.00019}
#'
#' Martinez, E. Z., Achcar, J. A. and Louzada-Neto, F. (2006). Estimators of sensitivity and specificity in the presence of verification bias: A Bayesian approach. *Computational Statistics & Data Analysis*, 51(2), 601–611. \doi{10.1016/j.csda.2005.12.021}
#'
#' Stan Development Team (2025). RStan: The R interface to Stan. R package. <https://CRAN.R-project.org/package=rstan>
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
                 if (isTRUE(sdata$J0 > 0L)) "q", if (sdata$n_beta > 0L) "beta",
                 # intercepts of verification, needed by compare_models()
                 if (isTRUE(sdata$mnar == 1L) && any(sdata$n1 + sdata$n0 > 0L)) "alpha")
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
  settings <- c(list(chains = chains, iter_warmup = iter_warmup,
                     iter_sampling = iter_sampling, seed = seed), list(...))
  dots <- list(...)
  # Saturated model: with few verified subjects of one class in a profile, the
  # posterior of logit(phi1) has a long tail against a steep side, which gives
  # divergent transitions at the default adapt_delta = 0.8
  if (S7::S7_inherits(model, saturated)) {
    if (backend == "rstan") {
      dots$control <- utils::modifyList(list(adapt_delta = 0.95), as.list(dots$control))
    } else if (is.null(dots$adapt_delta)) {
      dots$adapt_delta <- 0.95
    }
  }

  if (backend == "rstan") {
    stan_fit <- do.call(rstan::sampling, c(list(
      stanmodels[[program]], data = sdata, pars = variables,
      chains = chains, iter = iter_warmup + iter_sampling, warmup = iter_warmup,
      seed = seed, refresh = 0), dots))
    draws <- posterior::as_draws_array(rstan::extract(stan_fit, pars = variables,
                                                      permuted = FALSE))
  } else {
    stan_fit <- do.call(cmdstan_program(program)$sample, c(list(
      data = sdata, chains = chains, iter_warmup = iter_warmup,
      iter_sampling = iter_sampling, seed = seed, refresh = 0,
      show_messages = FALSE), dots))
    draws <- stan_fit$draws(variables)
  }

  baysesp_fit(model = model, data = data, draws = draws, backend = backend,
              stan_fit = stan_fit, settings = settings)
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
    stan_fit = S7::class_any,
    settings = S7::class_list
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
  d <- convergence(fit)
  sprintf("Convergence: max R-hat %s, min bulk ESS %s, min tail ESS %s, %s divergent transitions",
          formatC(d$max_rhat, digits = 3, format = "f"),
          formatC(d$min_ess_bulk, digits = 0, format = "f"),
          formatC(d$min_ess_tail, digits = 0, format = "f"), d$divergent)
}

# Worst convergence diagnostics as a one-row data frame (NA if not available).
convergence <- function(fit) {
  d <- posterior::summarise_draws(fit@draws, "rhat", "ess_bulk", "ess_tail")
  worst <- function(v, fun) {
    v <- v[is.finite(v)]
    if (length(v) == 0L) NA_real_ else fun(v)
  }
  data.frame(max_rhat = worst(d$rhat, max), min_ess_bulk = worst(d$ess_bulk, min),
             min_ess_tail = worst(d$ess_tail, min), divergent = divergent_transitions(fit))
}

#' Fit a different model to the same data
#'
#' Fits a new model to the data of a model already fitted with [baysesp()],
#' with the same backend and the same sampling settings (number of chains,
#' warmup and sampling iterations, seed, and the further arguments given to
#' `baysesp()`). Only the model changes: its class, its priors or its
#' verification mechanism. This is the usual way to fit the alternative models
#' given to [compare_models()].
#'
#' @param object A model fitted with [baysesp()].
#' @param model The new model: an object created by [saturated()],
#'   [loglinear()] or [random_effects()], with its priors and its verification
#'   mechanism (`verification = mar()` or `verification = mnar()`). Default
#'   `object@model`, the model of `object`, which fits the same model again.
#' @param ... Sampling settings that replace those of `object`, with the names
#'   of the arguments of [baysesp()]: `backend`, `chains`, `iter_warmup`,
#'   `iter_sampling`, `seed`, or arguments for the backend such as `cores`.
#'   Settings that are not given are those used for `object`; the seed, if
#'   `object` was fitted with `seed = NULL`, is the one drawn at that time.
#' @return A new `baysesp_fit` object, as returned by [baysesp()].
#' @seealso [baysesp()], [compare_models()]
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- simulate_study(n = 1000, prevalence = 0.2,
#'                     sensitivity = c(A = 0.9, B = 0.8),
#'                     specificity = c(A = 0.85, B = 0.95),
#'                     verification_prob = c(0.2, 0.9, 0.9, 1))
#' d <- baysesp_data(d, tests = c("A", "B"), verified = "verified", disease = "disease")
#' fit <- baysesp(d, seed = 1)
#' fit_ci <- update(fit, model = loglinear(interaction_order = 1))
#' fit_ci
#' }
#' @name update.baysesp_fit
#' @usage \method{update}{baysesp_fit}(object, model = object@model, ...)
NULL

S7::method(update_generic, baysesp_fit) <- function(object, model = object@model, ...) {
  args <- utils::modifyList(c(list(backend = object@backend), object@settings), list(...))
  do.call(baysesp, c(list(data = object@data, model = model), args))
}

#' Summary of a fitted model
#'
#' Collects in one object what describes a model fitted with [baysesp()]: the
#' model and its verification mechanism, the priors, the profiles of test
#' results with no verified subject, the convergence diagnostics of the
#' sampler and the posterior summaries of the accuracy of the tests.
#'
#' @details The convergence diagnostics are computed over all the variables
#'   saved in the draws (variables that are constant, such as the verified
#'   fraction of a never-verified profile, are skipped): the largest
#'   split-\eqn{\hat{R}}, which should be below 1.01, the smallest bulk and tail
#'   effective sample sizes, which should be at least 100 per chain (Vehtari et
#'   al., 2021), and the number of divergent transitions after warmup, which
#'   should be 0: divergences mean that the sampler could not explore some
#'   region of the posterior, and the draws may be biased.
#'
#' @param object A model fitted with [baysesp()].
#' @param level Probability of the credible intervals in `accuracy`, a single
#'   number in (0, 1). Default `0.95`.
#' @param ... Not used.
#' @return A list of class `fit_summary`, with a `print()` method, with the
#'   elements:
#'   * `model`: a description of the model, such as
#'     `"log-linear model of order 1"`.
#'   * `verification`: `"MAR"` or `"MNAR"`, followed by
#'     `"(one odds ratio per profile)"` with `by_profile = TRUE`.
#'   * `priors`: a data frame with the name of each argument that holds a
#'     prior (`argument`) and the prior (`prior`). The prior on the
#'     proportion of diseased in the never-verified profiles has one row per
#'     never-verified profile, named after it, such as
#'     `"unverified_prevalence_prior (0,1)"`. The priors of the
#'     verification mechanism are included when they are used, with the name
#'     of the argument of `mnar()`: `"mnar(odds_ratio_prior)"` with
#'     `verification = mnar()`, and `"mnar(intercepts_prior)"` with
#'     `verification = mnar()` in the log-linear and random effects models.
#'   * `never_verified`: the never-verified profiles, written as `"(0,1)"`,
#'     or a character vector of length 0.
#'   * `diagnostics`: a one-row data frame with `chains`, `draws` (per chain),
#'     `backend`, `max_rhat`, `min_ess_bulk`, `min_ess_tail` and `divergent`
#'     (number of divergent transitions; `NA` if not available).
#'   * `accuracy`: the result of [test_accuracy()], with the odds ratio of
#'     verification added under `verification = mnar()`.
#'   * `level`: the value of `level`.
#' @references
#' Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Bürkner, P.-C. (2021). Rank-normalization, folding, and localization: An improved \eqn{\hat{R}} for assessing convergence of MCMC (with discussion). *Bayesian Analysis*, 16(2), 667–718. \doi{10.1214/20-BA1221}
#' @seealso [baysesp()], [test_accuracy()]
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- simulate_study(n = 1000, prevalence = 0.2,
#'                     sensitivity = c(A = 0.9, B = 0.8),
#'                     specificity = c(A = 0.85, B = 0.95),
#'                     verification_prob = c(0.2, 0.9, 0.9, 1))
#' d <- baysesp_data(d, tests = c("A", "B"), verified = "verified", disease = "disease")
#' fit <- baysesp(d, seed = 1)
#' summary(fit)
#' }
#' @name summary.baysesp_fit
#' @usage \method{summary}{baysesp_fit}(object, level = 0.95, ...)
NULL

S7::method(summary_generic, baysesp_fit) <- function(object, level = 0.95, ...) {
  check_level(level)
  model <- object@model
  verification <- model@verification
  props <- S7::props(model)
  never <- object@data@never_verified
  priors <- unlist(lapply(names(props)[endsWith(names(props), "_prior")], function(p) {
    if (p != "unverified_prevalence_prior" || !any(never)) return(as_prior_list(props[[p]], p))
    # one row per never-verified profile, named after it
    q <- q_priors(unname(as_prior_list(props[[p]], p)), sum(never), p)
    stats::setNames(q, paste0(p, " ", profile_labels(object@data)[never]))
  }))
  if (S7::S7_inherits(verification, mnar)) {
    # named as given to mnar(), so they are not confused with the model's priors
    priors <- c(priors, list("mnar(odds_ratio_prior)" = verification@odds_ratio_prior),
                if (!S7::S7_inherits(model, saturated))
                  list("mnar(intercepts_prior)" = verification@intercepts_prior))
  }
  structure(list(
    model = sub("^<[^>]+> ", "", utils::capture.output(print(model))[1]),
    verification = paste0(toupper(class_label(verification)),
                          if (isTRUE(S7::prop_exists(verification, "by_profile")) &&
                              verification@by_profile) " (one odds ratio per profile)"),
    priors = data.frame(argument = names(priors),
                        prior = vapply(priors, prior_label, character(1)), row.names = NULL),
    never_verified = profile_labels(object@data)[object@data@never_verified],
    diagnostics = cbind(data.frame(chains = posterior::nchains(object@draws),
                                   draws = posterior::niterations(object@draws),
                                   backend = object@backend), convergence(object)),
    accuracy = rbind(test_accuracy(object, level), odds_ratio_summary(object, level)),
    level = level
  ), class = "fit_summary")
}

class_label <- function(x) sub("^baysesp::", "", class(x)[1])

#' @export
print.fit_summary <- function(x, ...) {
  d <- x$diagnostics
  cat("Model:          ", x$model, "\n",
      "Verification:   ", x$verification, "\n",
      "Never verified: ", if (length(x$never_verified) == 0L) "none" else
        paste(x$never_verified, collapse = " "), "\n\nPriors:\n", sep = "")
  for (i in seq_len(nrow(x$priors))) {
    cat("  ", format(x$priors$argument, width = max(nchar(x$priors$argument)))[i],
        "  ", x$priors$prior[i], "\n", sep = "")
  }
  cat("\nSampling: ", d$chains, " chains x ", d$draws, " draws (", d$backend, ")\n",
      "  max R-hat ", formatC(d$max_rhat, digits = 3, format = "f"),
      ", min bulk ESS ", formatC(d$min_ess_bulk, digits = 0, format = "f"),
      ", min tail ESS ", formatC(d$min_ess_tail, digits = 0, format = "f"),
      ", divergent transitions ", d$divergent, "\n\n",
      "Accuracy (", 100 * x$level, "% credible intervals):\n", sep = "")
  print(x$accuracy, row.names = FALSE, digits = 3)
  invisible(x)
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
