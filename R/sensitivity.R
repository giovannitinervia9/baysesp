#' @include accuracy.R stan_data.R
NULL

# Sensitivity analysis for the saturated model ----------------------------------
#
# In the saturated model the parameters that the data cannot estimate, q (the
# proportion of diseased in never-verified profiles) and beta (MNAR), do not
# enter the likelihood and are independent a priori of the identified ones
# (omega, v, phi1). Their posterior is their prior, independent of the
# posterior of the identified parameters. So the analysis under different
# assumptions needs no new fit: keep the draws of omega, v and phi1, replace
# those of q and beta, and recompute (specification, sections 7.2 and 10.2):
#   psi_t = v_t phi1_t + (1 - v_t) expit(logit(phi1_t) - beta)  (verified profiles)
#   psi_t = q_t                                                 (never verified)
#   prevalence = sum(omega * psi), theta_1 = omega * psi / prevalence, ...

#' Sensitivity analysis for the assumptions the data cannot check
#'
#' Shows how the accuracy of the tests changes with the two assumptions of the
#' saturated model that the data cannot check: the odds ratio of verification
#' between diseased and healthy subjects (see [mnar()]) and the prior on the
#' proportion of diseased in the never-verified profiles. For the saturated
#' model no new fit is needed: these parameters do not enter the likelihood,
#' so the draws of the identified parameters are reused and only the
#' assumptions are changed.
#'
#' For the structured models ([loglinear()]) the odds ratio enters the
#' likelihood, so the model is fitted again for each value, with the odds ratio
#' fixed; the never-verified profiles have no separate prior to change.
#'
#' The odds ratio is fixed at each value of `odds_ratios`: with a grid of
#' values this is a tipping-point analysis, which shows at which odds ratio a
#' conclusion changes. With `odds_ratios = NULL` the odds ratio of the fitted
#' model is kept (1 for [mar()], the prior of [mnar()] otherwise). The new
#' priors in `unverified_priors` are sampled with the random number generator
#' of R: set the seed for reproducible results.
#'
#' @param fit A model fitted with [baysesp()].
#' @param odds_ratios Values of the odds ratio of verification between diseased
#'   and healthy subjects, \eqn{\exp(\beta)}; 1 is MAR. `NULL` keeps the
#'   assumption of the fitted model.
#' @param unverified_priors List of alternative priors on the proportion of
#'   diseased in the never-verified profiles. Each element is a prior, or a list
#'   with one prior per never-verified profile. `NULL` keeps the prior of the
#'   fitted model.
#' @param rules Rules, see [rule_accuracy()]. By default the individual tests,
#'   as in [test_accuracy()].
#' @inheritParams test_accuracy
#' @param ... For the structured models, arguments passed to [baysesp()] for
#'   the new fits, such as `chains` and `seed`.
#' @return A data frame with the columns `odds_ratio` (`NA` when it follows the
#'   prior of an [mnar()] model) and `unverified_prior`, which describe the
#'   assumptions, followed by the columns of [test_accuracy()] or
#'   [rule_accuracy()], without `identified`.
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- simulate_study(n = 2000, prevalence = 0.2,
#'                     sensitivity = c(A = 0.9, B = 0.8),
#'                     specificity = c(A = 0.85, B = 0.95),
#'                     verification_prob = c(0.2, 0.9, 0.9, 1), odds_ratio = 2)
#' d <- baysesp_data(d, tests = c("A", "B"), verified = "verified", disease = "disease")
#' fit <- baysesp(d, seed = 1)
#' sensitivity_analysis(fit, odds_ratios = c(1, 1.5, 2, 3))
#' }
#' @export
sensitivity_analysis <- function(fit, odds_ratios = NULL, unverified_priors = NULL,
                                 rules = NULL, level = 0.95, ...) {
  if (!S7::S7_inherits(fit, baysesp_fit)) {
    stop("`fit` must be a model fitted with baysesp().", call. = FALSE)
  }
  check_level(level)
  if (!is.null(odds_ratios) && (!is.numeric(odds_ratios) || anyNA(odds_ratios) ||
                                any(odds_ratios <= 0) || any(!is.finite(odds_ratios)))) {
    stop("`odds_ratios` must contain positive finite numbers.", call. = FALSE)
  }
  if (S7::S7_inherits(fit@model, loglinear)) {
    return(refit_sensitivity(fit, odds_ratios, unverified_priors, rules, level, ...))
  }
  check_saturated_fit(fit, "sensitivity_analysis()")
  t <- transparent_draws(fit)
  S <- nrow(t$omega)
  never <- fit@data@never_verified
  J0 <- sum(never)

  if (is.null(odds_ratios)) {
    betas <- list(t$beta)
    or_labels <- if (S7::S7_inherits(fit@model@verification, mar)) 1 else NA_real_
  } else {
    betas <- lapply(log(odds_ratios), function(b) matrix(b, S, ncol(t$omega)))
    or_labels <- odds_ratios
  }

  if (is.null(unverified_priors)) {
    qs <- list(t$q)
    q_labels <- if (J0 == 0L) NA_character_ else
      prior_labels(as_prior_list(fit@model@unverified_prevalence_prior, ""))
  } else {
    if (J0 == 0L) {
      stop("`unverified_priors` is given, but every profile has verified subjects.",
           call. = FALSE)
    }
    if (S7::S7_inherits(unverified_priors, prior_distribution)) {
      unverified_priors <- list(unverified_priors)
    }
    if (!is.list(unverified_priors)) {
      stop("`unverified_priors` must be a list of priors.", call. = FALSE)
    }
    priors <- lapply(unverified_priors, function(p) {
      p <- as_prior_list(p, "unverified_priors")
      for (name in names(p)) check_prior(p[[name]], "probability", name)
      q_priors(unname(p), J0, "unverified_priors")
    })
    qs <- lapply(priors, function(p) vapply(p, prior_sample, numeric(S), n = S))
    q_labels <- vapply(priors, prior_labels, character(1))
  }

  grid <- expand.grid(q = seq_along(qs), beta = seq_along(betas))
  out <- lapply(seq_len(nrow(grid)), function(i) {
    scenario <- fit
    scenario@draws <- core_from_transparent(t, betas[[grid$beta[i]]], qs[[grid$q[i]]], never)
    acc <- if (is.null(rules)) test_accuracy(scenario, level) else
      rule_accuracy(scenario, rules, level)
    acc$identified <- NULL
    cbind(odds_ratio = or_labels[grid$beta[i]], unverified_prior = q_labels[grid$q[i]], acc)
  })
  do.call(rbind, out)
}

# Structured models: one new fit per odds ratio, with beta fixed in Stan.
refit_sensitivity <- function(fit, odds_ratios, unverified_priors, rules, level, ...) {
  if (!is.null(unverified_priors)) {
    stop("In the structured models the never-verified profiles have no separate ",
         "prior: `unverified_priors` applies only to saturated().", call. = FALSE)
  }
  if (is.null(odds_ratios)) {
    stop("Give the values of the odds ratio with `odds_ratios`.", call. = FALSE)
  }
  model <- fit@model
  if (!S7::S7_inherits(model@verification, mnar)) model@verification <- mnar()
  args <- utils::modifyList(list(chains = 4, iter_warmup = 1000, iter_sampling = 1000,
                                 seed = NULL), list(...))
  out <- lapply(odds_ratios, function(or) {
    sdata <- stan_data(model, model@verification, fit@data)
    sdata$n_beta <- 0L
    sdata$beta_value <- log(or)
    refit <- do.call(fit_stan, c(list(model, fit@data, sdata, fit@backend), args))
    acc <- if (is.null(rules)) test_accuracy(refit, level) else rule_accuracy(refit, rules, level)
    acc$identified <- NULL
    cbind(odds_ratio = or, unverified_prior = NA_character_, acc)
  })
  do.call(rbind, out)
}

check_saturated_fit <- function(fit, what) {
  if (!S7::S7_inherits(fit, baysesp_fit)) {
    stop("`fit` must be a model fitted with baysesp().", call. = FALSE)
  }
  if (!S7::S7_inherits(fit@model, saturated)) {
    stop(what, " is available only for the saturated model.", call. = FALSE)
  }
}

prior_labels <- function(priors) {
  paste(vapply(priors, prior_label, character(1)), collapse = "; ")
}

# Draws of the transparent parametrization as draws x profiles matrices; q is
# draws x never-verified profiles, beta draws x profiles (0 under MAR).
transparent_draws <- function(fit) {
  m <- unclass(posterior::as_draws_matrix(fit@draws))
  J <- nrow(fit@data@profiles)
  never <- fit@data@never_verified
  get <- function(name, n) {
    if (n == 0L) return(matrix(0, nrow(m), 0L))
    unname(m[, paste0(name, "[", seq_len(n), "]"), drop = FALSE])
  }
  beta <- matrix(0, nrow(m), J)
  n_beta <- sum(startsWith(colnames(m), "beta["))
  if (n_beta > 0L) beta[, !never] <- get("beta", n_beta)
  list(omega = get("omega", J), v = get("v", J), phi1 = get("phi1", J),
       q = get("q", sum(never)), beta = beta)
}

# Proportion of diseased in each profile, psi, for given beta (draws x profiles);
# 0 in the never-verified profiles.
profile_prevalence <- function(t, beta, never) {
  psi <- t$v * t$phi1 + (1 - t$v) * stats::plogis(stats::qlogis(t$phi1) - beta)
  psi[, never] <- 0
  psi
}

# Draws of prevalence and profile distributions for given beta and q.
core_from_transparent <- function(t, beta, q, never) {
  psi <- profile_prevalence(t, beta, never)
  psi[, never] <- q
  w <- t$omega * psi
  prevalence <- rowSums(w)
  J <- ncol(w)
  x <- cbind(prevalence, (t$omega - w) / (1 - prevalence), w / prevalence)
  colnames(x) <- c("prevalence", paste0("theta_healthy[", seq_len(J), "]"),
                   paste0("theta_diseased[", seq_len(J), "]"))
  posterior::as_draws_matrix(x)
}

# Draws from a prior on a probability; priors on the real line act on the logit.
prior_sample <- function(prior, n) {
  switch(prior@family,
    beta = stats::rbeta(n, prior@shape1, prior@shape2),
    normal = stats::plogis(stats::rnorm(n, prior@mean, prior@sd)),
    student = stats::plogis(prior@location + prior@scale * stats::rt(n, prior@df)),
    stop("Cannot sample from a ", prior@family, " prior.", call. = FALSE)
  )
}
