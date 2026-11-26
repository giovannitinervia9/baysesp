#' @include fit.R sensitivity.R plots.R generics.R
NULL

# Comparison of models fitted to the same data ------------------------------------
#
# Every model gives a probability to the observed cells: for each profile t,
# verified diseased, verified healthy and unverified (the verified cells of
# never-verified profiles are structurally empty and left out). The likelihood
# includes the verification part also under MAR, with a free verification
# probability per profile, so that MAR and MNAR models can be compared.
# The observation is the subject, and the subjects of a cell are identical:
# leave-one-out has at most 3J distinct terms, each counted as many times as
# the subjects in its cell. The saturated model has an exact leave-one-out in
# closed form; the structured models use PSIS-LOO, with an exact refit for the
# cells with Pareto k > 0.7.

#' Compare models fitted to the same data
#'
#' Compares models fitted with [baysesp()] to the same data by their expected
#' log predictive density for a new subject, estimated by leave-one-out
#' cross-validation (LOO), with the widely applicable information criterion
#' (WAIC) as a second estimate, and checks the fit of each model with posterior
#' predictive checks on the counts of the observed table.
#'
#' @details
#' **What is compared.** The data are the counts of the cells of the observed
#' table: for each response profile \eqn{t} of the tests, the verified
#' diseased \eqn{n_{t,1}}, the verified healthy \eqn{n_{t,0}} and the
#' unverified \eqn{u_t} (the verified cells of the profiles with no verified
#' subject are empty by design and are left out). Every model gives each cell
#' a probability: with prevalence \eqn{\pi}, distributions of the profiles
#' \eqn{\theta_{d,t}} among the diseased (\eqn{d = 1}) and the healthy
#' (\eqn{d = 0}) and probabilities of verification
#' \eqn{\lambda_{t,d} = P(V = 1 \mid T = t, D = d)},
#' \deqn{P(t, V = 1, D = d) = \pi_d\, \theta_{d,t}\, \lambda_{t,d}, \qquad
#' P(t, V = 0) = \sum_d \pi_d\, \theta_{d,t}\, (1 - \lambda_{t,d}),}
#' with \eqn{\pi_1 = \pi} and \eqn{\pi_0 = 1 - \pi}.
#'
#' **The part of verification.** The probability of the observed table
#' includes which subjects were verified, also under `verification = mar()`,
#' where this part is ignored when the model is fitted: this makes the models
#' with `verification = mar()` and `verification = mnar()` comparable. With
#' MAR, \eqn{\lambda_{t,1} = \lambda_{t,0} = v_t}, a free probability of
#' verification for each profile with verified subjects, independent of the
#' other parameters; the same \eqn{v_t} is the verified fraction of the
#' saturated model, also with `verification = mnar()`, where the odds ratio
#' does not enter the likelihood. In all these models, \eqn{v_t} is given the
#' prior `verified_prior`, the same for every model, in place of the
#' `verified_prior` of a saturated model. Its contribution is then the same in
#' every model and cancels in the differences between them: these models are
#' compared only on how they describe the tests. Its posterior is
#' \eqn{p(v_t) \propto \textrm{prior}(v_t)\, v_t^{n_{t,1} + n_{t,0}} (1 - v_t)^{u_t}},
#' a Beta distribution for a Beta prior; its draws are taken in R (set the
#' seed for reproducible results), and its leave-one-out predictions are exact.
#' Under `verification = mnar()` in the log-linear and random effects models,
#' verification depends on the disease and cannot be separated: it follows the
#' priors of the model, `intercepts_prior` and `odds_ratio_prior` of `mnar()`.
#' To compare these models with the others, choose a `verified_prior` close to
#' the prior that `intercepts_prior` gives to verification: the default
#' `verified_prior`, uniform, is close to the default
#' `intercepts_prior = prior_normal(0, 1.5)`, and `verified_prior =
#' prior_normal(0, 2.5)` matches `intercepts_prior = prior_normal(0, 2.5)`.
#' Then the comparison is not driven by different priors on verification.
#'
#' **Leave-one-out.** The unit left out is the subject, not the cell. The
#' subjects of the same cell are identical, so there are at most as many
#' distinct leave-one-out predictions as cells, and the contribution of a cell
#' is multiplied by its count. The expected log predictive density is
#' \deqn{\mathrm{elpd}_{\mathrm{loo}} = \sum_c n_c \log p(y_c \mid y_{-c}),}
#' where \eqn{n_c} is the count of cell \eqn{c} and
#' \eqn{p(y_c \mid y_{-c})} the posterior predictive probability of a subject
#' in cell \eqn{c} given all the other subjects (Vehtari, Gelman and Gabry,
#' 2017). Its standard error is \eqn{\sqrt{N}} times the standard deviation of
#' the \eqn{N} contributions of the subjects. Higher values mean better
#' predictions; a difference smaller than about twice its standard error is
#' not clear evidence.
#'
#' For the saturated model the leave-one-out predictions are exact and in
#' closed form: the posterior of the profile frequencies is Dirichlet and the
#' posteriors of \eqn{v_t} and of the proportion of diseased among the
#' verified \eqn{\phi^{(1)}_t} are independent, so the predictions are ratios
#' of Dirichlet-multinomial and Beta-binomial marginal likelihoods (with
#' `verified_prior` for \eqn{v_t}). With a
#' prior on the logit scale (normal or Student-t) the one-dimensional marginal
#' likelihoods are computed by numerical integration. The proportion of
#' diseased in the never-verified profiles and the odds ratio of verification
#' do not enter the likelihood of the saturated model and do not change its
#' predictions.
#'
#' For the log-linear and random effects models, leave-one-out is estimated by
#' Pareto smoothed importance sampling (PSIS-LOO) from the posterior draws
#' (Vehtari et al., 2024). The estimated shape \eqn{k} of the Pareto
#' distribution says whether the estimate of a cell is reliable: above 0.7 it
#' is not, and with `refit = TRUE` the model is fitted again without one
#' subject of that cell, with the sampling settings of the original fit, and
#' the prediction is computed exactly from the new draws. There are at most as
#' many new fits as cells. Small cells, such as the verified diseased among
#' the subjects with all tests negative, are the influential ones.
#'
#' **WAIC.** The widely applicable information criterion (Watanabe, 2010) is
#' an asymptotic approximation of leave-one-out. It is computed cell by cell,
#' as the log of the posterior mean of the probability of the cell minus its
#' posterior variance on the log scale, and multiplied by the counts. When it
#' differs from the leave-one-out estimate by more than 1, a warning is given:
#' the approximation is not reliable for these data, usually because of small
#' cells, and the leave-one-out estimate should be used.
#'
#' **Posterior predictive checks.** For each posterior draw a table with the
#' same number of subjects is drawn from the probabilities of the cells
#' (Gelman, Meng and Stern, 1996). Each observed count is compared with the
#' distribution of its replicates, with a two-sided posterior predictive
#' p-value, \eqn{2 \min\{P(y^{\mathrm{rep}}_c \ge y_c), P(y^{\mathrm{rep}}_c \le y_c)\}}
#' (at most 1). The whole table is compared with the chi-square discrepancy
#' \eqn{X^2 = \sum_c (y_c - E_c)^2 / E_c}, with \eqn{E_c} the expected count
#' under the draw: the p-value is the proportion of draws in which the
#' discrepancy of the replicated table is at least that of the observed table.
#' Small p-values mean that the model does not reproduce the data.
#'
#' **What the data cannot tell.** All the models are compared on the observed
#' cells only. Two saturated models that differ only in the prior on the
#' proportion of diseased in the never-verified profiles
#' (`unverified_prevalence_prior`), in the verification mechanism
#' (`verification = mar()` or `verification = mnar()`) or in their own
#' `verified_prior`, replaced here by the common one, give the same
#' probability to the observed data: the comparison cannot tell them apart, and
#' a warning says so. The effect of these choices is shown by varying them,
#' as `sensitivity_analysis()` does. In the same way, a structured model that
#' predicts the observed cells well is not confirmed in what it assumes about
#' the unobserved disease status of the unverified subjects.
#'
#' @param ... Two or more models fitted with [baysesp()] to the same data. A
#'   single model is also accepted: then only its leave-one-out estimate and
#'   its checks are computed. The names of the arguments, if given (such as
#'   `compare_models(saturated = fit, independence = fit_ci)`), are the labels
#'   of the models; otherwise the expressions are used. Models fitted to
#'   different data give an error.
#' @param refit `TRUE` or `FALSE`, a single value. With `TRUE` (default), for
#'   the log-linear and random effects models, each cell with Pareto
#'   \eqn{k > 0.7} is estimated exactly with a new fit; with `FALSE`, no new
#'   fit is run and a warning lists the models with unreliable cells. A new
#'   fit is also skipped, with a warning, when with `verification = mnar()`
#'   leaving out the subject would leave its profile with no verified subject.
#' @param verified_prior Prior on the probability of verification of each
#'   profile, used for the part of verification of every model where it can be
#'   separated: the saturated models and the log-linear and random effects
#'   models with `verification = mar()` (see Details). A single prior on a
#'   probability: `prior_beta()`, or `prior_normal()` or `prior_student()`,
#'   which act on its logit. Default `prior_beta(0.5, 2)`, the uniform
#'   distribution. It replaces the `verified_prior` of the saturated models
#'   only within the comparison; their estimates do not change.
#' @param level Probability of the intervals of the replicated counts in
#'   `checks`, a single number in (0, 1). Default `0.95`.
#' @return A list of class `model_comparison`, with `print()` and `plot()`
#'   methods and two data frames:
#'   * `table`: one row per model, ordered by decreasing `elpd_loo`, with
#'     `model` (the label), `elpd_loo` and its standard error `se_loo`,
#'     `elpd_diff` (difference from the first model, 0 for the first) and its
#'     standard error `se_diff`, `p_loo` (effective number of parameters, the
#'     log predictive density of the data minus `elpd_loo`), `elpd_waic` and
#'     `p_waic`, `loo_method` (`"exact"`, `"psis"` or `"psis"` followed by the
#'     number of new fits), `high_k` (number of cells with Pareto
#'     \eqn{k > 0.7}; `NA` for the saturated model, whose leave-one-out is
#'     exact and has no Pareto \eqn{k}, printed as `-`) and `ppc_pvalue` (the
#'     p-value of the chi-square discrepancy).
#'   * `checks`: one row per model and cell, with `model`, `profile` (written
#'     as `"(0,1)"`), `cell` (`"diseased"`, `"healthy"` or `"unverified"`),
#'     `observed` (the count), `expected` (posterior predictive mean of the
#'     count), `lower` and `upper` (interval of the replicated count with
#'     probability `level`) and `p_value` (two-sided posterior predictive
#'     p-value of the cell).
#' @references
#' Gelman, A., Meng, X.-L. and Stern, H. (1996). Posterior predictive assessment of model fitness via realized discrepancies. *Statistica Sinica*, 6(4), 733–807. <https://www3.stat.sinica.edu.tw/statistica/j6n4/j6n41/j6n41.htm>
#'
#' Vehtari, A., Gelman, A. and Gabry, J. (2017). Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. *Statistics and Computing*, 27(5), 1413–1432. \doi{10.1007/s11222-016-9696-4}
#'
#' Vehtari, A., Simpson, D., Gelman, A., Yao, Y. and Gabry, J. (2024). Pareto smoothed importance sampling. *Journal of Machine Learning Research*, 25(72), 1–58. <https://jmlr.org/papers/v25/19-556.html>
#'
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory. *Journal of Machine Learning Research*, 11(116), 3571–3594. <https://jmlr.org/papers/v11/watanabe10a.html>
#' @seealso [baysesp()], [sensitivity_analysis()]
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
#' cmp <- compare_models(saturated = fit, independence = fit_ci)
#' cmp
#' plot(cmp)
#' }
#' @export
compare_models <- function(..., refit = TRUE, verified_prior = prior_beta(0.5, 2),
                           level = 0.95) {
  fits <- list(...)
  if (length(fits) == 0L || !all(vapply(fits, S7::S7_inherits, logical(1), baysesp_fit))) {
    stop("Give one or more models fitted with baysesp().", call. = FALSE)
  }
  if (!isTRUE(refit) && !isFALSE(refit)) {
    stop("`refit` must be TRUE or FALSE.", call. = FALSE)
  }
  check_prior(verified_prior, "probability", "verified_prior")
  check_level(level)
  labels <- vapply(as.list(substitute(list(...)))[-1], deparse1, character(1))
  if (!is.null(names(fits))) labels[names(fits) != ""] <- names(fits)[names(fits) != ""]
  labels <- make.unique(labels)
  data <- as.data.frame(fits[[1]]@data)
  if (!all(vapply(fits, function(f) identical(as.data.frame(f@data), data), logical(1)))) {
    stop("The models must be fitted to the same data.", call. = FALSE)
  }
  same_likelihood_warning(fits, labels)

  cells <- model_cells(fits[[1]]@data)
  results <- Map(function(f, label) {
    evaluate_model(f, label, cells, refit, verified_prior, level)
  }, fits, labels)
  counts <- cells$count[cells$count > 0L]
  elpd <- vapply(results, function(r) sum(counts * r$elpd), numeric(1))
  best <- which.max(elpd)
  table <- do.call(rbind, lapply(results, function(r) {
    loo <- weighted_total(r$elpd, counts)
    diff <- weighted_total(r$elpd - results[[best]]$elpd, counts)
    waic <- sum(counts * r$waic)
    data.frame(model = r$label, elpd_loo = loo[1], se_loo = loo[2],
               elpd_diff = diff[1], se_diff = diff[2],
               p_loo = sum(counts * r$lpd) - loo[1],
               elpd_waic = waic, p_waic = sum(counts * r$p_waic),
               loo_method = r$method, high_k = r$high_k, ppc_pvalue = r$ppc_pvalue)
  }))
  table <- table[order(-table$elpd_loo), ]
  rownames(table) <- NULL

  gap <- abs(table$elpd_waic - table$elpd_loo) > 1
  if (any(gap)) {
    warning("WAIC differs from LOO by more than 1 for ",
            paste(table$model[gap], collapse = ", "),
            ": use LOO, WAIC is not reliable for these data.", call. = FALSE)
  }
  checks <- do.call(rbind, lapply(results, `[[`, "checks"))
  rownames(checks) <- NULL
  structure(list(table = table, checks = checks), class = "model_comparison")
}

# Cells of the observed table: profile, type and count, without the verified
# cells of the never-verified profiles.
model_cells <- function(data) {
  J <- length(data@n_diseased)
  cells <- data.frame(profile = rep(seq_len(J), each = 3L),
                      cell = rep(c("diseased", "healthy", "unverified"), J),
                      count = c(rbind(data@n_diseased, data@n_healthy, data@n_unverified)))
  never <- data@never_verified[cells$profile]
  cells <- cells[!never | cells$cell == "unverified", ]
  cells$label <- profile_labels(data)[cells$profile]
  rownames(cells) <- NULL
  cells
}

# Probability of each cell in each draw (draws x cells). Where verification
# can be separated (saturated models, structured models under MAR), v_t is
# drawn from its posterior under the common `verified_prior`. `never`:
# profiles with no verified subject in the original data, whose verified
# cells have probability 0 (a refit may have fewer verified profiles).
cell_probs <- function(fit, cells, verified_prior, never = fit@data@never_verified) {
  d <- fit@data
  J <- length(never)
  if (S7::S7_inherits(fit@model, saturated)) {
    t <- transparent_draws(fit)
    v <- verified_draws(verified_prior, d, nrow(t$omega), never)
    p1 <- t$omega * v * t$phi1
    p0 <- t$omega * v * (1 - t$phi1)
    pu <- t$omega * (1 - v)
  } else {
    core <- core_draws(fit)
    j1 <- core$prevalence * core$theta_diseased
    j0 <- (1 - core$prevalence) * core$theta_healthy
    S <- nrow(j1)
    if (S7::S7_inherits(fit@model@verification, mnar)) {
      m <- unclass(posterior::as_draws_matrix(fit@draws))
      verified <- !d@never_verified
      alpha <- beta <- matrix(0, S, J)
      alpha[, verified] <- m[, paste0("alpha[", seq_len(sum(verified)), "]")]
      b <- m[, grep("^beta\\[", colnames(m)), drop = FALSE]
      if (ncol(b) == 1L) beta[] <- b[, 1] else beta[, verified] <- b
      l1 <- stats::plogis(alpha + beta)
      l0 <- stats::plogis(alpha)
    } else {
      l1 <- l0 <- verified_draws(verified_prior, d, S, never)
    }
    l1[, never] <- 0
    l0[, never] <- 0
    p1 <- j1 * l1
    p0 <- j0 * l0
    pu <- j1 * (1 - l1) + j0 * (1 - l0)
  }
  type <- match(cells$cell, c("diseased", "healthy", "unverified"))
  cbind(p1, p0, pu)[, (type - 1L) * J + cells$profile, drop = FALSE]
}

# Posterior draws (draws x profiles) of the free probability of verification
# of each profile, with `prior` and the counts of `d`; 0 in the profiles `never`.
verified_draws <- function(prior, d, S, never) {
  s <- d@n_diseased + d@n_healthy
  v <- vapply(seq_along(s), function(j) {
    if (never[j]) numeric(S) else posterior_sample(prior, s[j], d@n_unverified[j], S)
  }, numeric(S))
  matrix(v, S)
}

# n draws of a probability with prior `prior` after s successes and f failures:
# Beta in closed form, otherwise from a fine grid on the logit scale.
posterior_sample <- function(prior, s, f, n) {
  if (prior@family == "beta") return(stats::rbeta(n, prior@shape1 + s, prior@shape2 + f))
  # ponytail: grid of +-12 sd around the mode; widen it if heavy-tailed priors
  # with almost no data ever matter here
  p <- logit_posterior(prior, s, f)
  x <- seq(p$mode - 12 * p$sd, p$mode + 12 * p$sd, length.out = 4001)
  stats::plogis(sample(x, n, replace = TRUE, prob = exp(p$log_density(x) - p$top)))
}

# LOO, WAIC and posterior predictive checks of one model.
evaluate_model <- function(fit, label, cells, refit, verified_prior, level) {
  probs <- cell_probs(fit, cells, verified_prior)
  observed <- cells$count > 0L
  ll <- log(probs[, observed, drop = FALSE])
  lpd <- log(colMeans(exp(ll)))
  p_waic <- apply(ll, 2, stats::var)
  high_k <- NA_integer_
  if (S7::S7_inherits(fit@model, saturated)) {
    elpd <- exact_loo(fit, cells[observed, ], verified_prior)
    method <- "exact"
  } else {
    psis <- psis_loo(ll, posterior::nchains(fit@draws))
    elpd <- psis$elpd
    high <- which(psis$k > 0.7)
    high_k <- length(high)
    method <- "psis"
    if (refit && length(high) > 0L) {
      message("Refitting ", label, " for ", length(high),
              " cells with Pareto k > 0.7.")
      new <- vapply(high, function(i) {
        refit_loo(fit, cells[observed, ][i, ], verified_prior)
      }, numeric(1))
      elpd[high[!is.na(new)]] <- new[!is.na(new)]
      if (any(!is.na(new))) method <- paste0("psis, ", sum(!is.na(new)), " refits")
      high <- high[is.na(new)]
    }
    if (length(high) > 0L) {
      warning("PSIS-LOO is not reliable for ", length(high), " cells of ", label,
              " (Pareto k > 0.7)", if (!refit) ": use refit = TRUE", ".", call. = FALSE)
    }
  }
  check <- predictive_check(probs, cells$count, level)
  list(label = label, elpd = elpd, lpd = lpd, waic = lpd - p_waic, p_waic = p_waic,
       method = method, high_k = high_k, ppc_pvalue = check$p_value,
       checks = data.frame(model = label, profile = cells$label, cell = cells$cell,
                           check$cells))
}

# PSIS-LOO of each column of a draws x cells log-likelihood matrix, with the
# draws ordered by chain.
psis_loo <- function(ll, chains) {
  chain <- rep(seq_len(chains), each = nrow(ll) / chains)
  r_eff <- loo::relative_eff(exp(ll), chain_id = chain)
  x <- suppressWarnings(loo::loo(ll, r_eff = r_eff))
  list(elpd = unname(x$pointwise[, "elpd_loo"]), k = unname(x$diagnostics$pareto_k))
}

# Exact leave-one-out of one subject in cell `cell`: fit again without it and
# average the probability of the cell over the new draws. NA (with a warning)
# when the refit would lose the intercept of verification of the profile.
refit_loo <- function(fit, cell, verified_prior) {
  d <- fit@data
  j <- cell$profile
  verified_cell <- cell$cell != "unverified"
  if (verified_cell && S7::S7_inherits(fit@model@verification, mnar) &&
      d@n_diseased[j] + d@n_healthy[j] == 1L) {
    warning("Cell ", cell$label, " ", cell$cell, ": no new fit, because without its ",
            "subject the profile has no verified subject.", call. = FALSE)
    return(NA_real_)
  }
  slot <- switch(cell$cell, diseased = "n_diseased", healthy = "n_healthy",
                 unverified = "n_unverified")
  counts <- S7::prop(d, slot)
  counts[j] <- counts[j] - 1L
  S7::prop(d, slot) <- counts
  new <- do.call(baysesp, c(list(data = d, model = fit@model, backend = fit@backend),
                            fit@settings))
  log(mean(cell_probs(new, cell, verified_prior, never = fit@data@never_verified)))
}

# Exact leave-one-out of the saturated model, one value per cell: Dirichlet-
# multinomial for the profile, Beta-binomial (or numerical integration) for
# verification, with `verified_prior`, and for the proportion of diseased
# among the verified.
exact_loo <- function(fit, cells, verified_prior) {
  d <- fit@data
  model <- fit@model
  n1 <- d@n_diseased
  n0 <- d@n_healthy
  u <- d@n_unverified
  N <- n1 + n0 + u
  a <- rep_len(model@frequency_prior@concentration, length(N))
  frequency <- log(a + N - 1) - log(sum(a) + sum(N) - 1)
  verified <- n1 + n0
  # log predictive of one success (or failure) given the others
  term <- function(prior, s, f, success) {
    log_marginal(prior, s, f) - log_marginal(prior, s - success, f - (1 - success))
  }
  vapply(seq_len(nrow(cells)), function(i) {
    j <- cells$profile[i]
    v <- verified_prior
    p <- model@profile_prevalence_prior
    frequency[j] + switch(cells$cell[i],
      diseased = term(v, verified[j], u[j], 1) + term(p, n1[j], n0[j], 1),
      healthy = term(v, verified[j], u[j], 1) + term(p, n1[j], n0[j], 0),
      unverified = if (verified[j] > 0L) term(v, verified[j], u[j], 0) else 0)
  }, numeric(1))
}

# log of the integral of p^s (1 - p)^f under a prior on a probability: closed
# form for a Beta prior, numerical integration on the logit scale otherwise.
log_marginal <- function(prior, s, f) {
  if (prior@family == "beta") {
    return(lbeta(prior@shape1 + s, prior@shape2 + f) - lbeta(prior@shape1, prior@shape2))
  }
  if (s + f == 0) return(0)
  p <- logit_posterior(prior, s, f)
  f_int <- function(x) exp(p$log_density(x) - p$top)
  # three pieces, so that a sharp peak is not missed
  cuts <- c(-Inf, p$mode - 10 * p$sd, p$mode + 10 * p$sd, Inf)
  total <- sum(vapply(1:3, function(i) {
    stats::integrate(f_int, cuts[i], cuts[i + 1L], rel.tol = 1e-10)$value
  }, numeric(1)))
  log(total) + p$top
}

# Unnormalised log posterior on the logit scale of a probability with a prior
# on its logit (normal or Student-t), after s successes and f failures, with
# its mode, its value at the mode and a scale from the curvature there.
logit_posterior <- function(prior, s, f) {
  log_prior <- switch(prior@family,
    normal = function(x) stats::dnorm(x, prior@mean, prior@sd, log = TRUE),
    student = function(x) stats::dt((x - prior@location) / prior@scale, prior@df,
                                    log = TRUE) - log(prior@scale))
  g <- function(x) s * stats::plogis(x, log.p = TRUE) +
    f * stats::plogis(-x, log.p = TRUE) + log_prior(x)
  mode <- stats::optimize(g, c(-50, 50), maximum = TRUE)$maximum
  top <- g(mode)
  h <- 1e-3
  sd <- 1 / sqrt(max(-(g(mode + h) - 2 * top + g(mode - h)) / h^2, 1e-8))
  list(log_density = g, mode = mode, top = top, sd = sd)
}

# Total over subjects of per-cell values `e`, with counts `n`, and its standard
# error: sqrt(N) times the standard deviation over the N subjects.
weighted_total <- function(e, n) {
  N <- sum(n)
  m <- sum(n * e) / N
  c(sum(n * e), sqrt(N * sum(n * (e - m)^2) / (N - 1)))
}

# Posterior predictive checks: replicated tables, one per draw.
predictive_check <- function(probs, y, level) {
  N <- sum(y)
  replicates <- vapply(seq_len(nrow(probs)), function(s) {
    as.numeric(stats::rmultinom(1L, N, probs[s, ]))
  }, numeric(length(y)))
  replicates <- matrix(replicates, nrow = length(y))
  expected <- t(probs) * N
  chisq <- function(x) colSums((x - expected)^2 / expected)
  tail <- (1 - level) / 2
  list(
    p_value = mean(chisq(replicates) >= chisq(matrix(y, length(y), nrow(probs)))),
    cells = data.frame(
      observed = y,
      expected = rowMeans(replicates),
      lower = apply(replicates, 1, stats::quantile, probs = tail, names = FALSE),
      upper = apply(replicates, 1, stats::quantile, probs = 1 - tail, names = FALSE),
      p_value = pmin(1, 2 * pmin(rowMeans(replicates >= y), rowMeans(replicates <= y)))
    )
  )
}

# Saturated models that differ only in what does not enter the likelihood.
same_likelihood_warning <- function(fits, labels) {
  key <- function(f) {
    m <- f@model
    if (!S7::S7_inherits(m, saturated)) return(NA_character_)
    # their own verified_prior is replaced by the common one
    paste(deparse(list(m@frequency_prior@concentration,
                       prior_label(m@profile_prevalence_prior))), collapse = "")
  }
  keys <- vapply(fits, key, character(1))
  for (k in unique(keys[!is.na(keys) & duplicated(keys)])) {
    warning("Models ", paste(labels[keys %in% k], collapse = ", "),
            " give the same probability to the observed data: they differ at most in ",
            "the prior on the never-verified profiles, in the verification mechanism ",
            "or in verified_prior, which do not change the comparison. The comparison ",
            "cannot tell them apart: use sensitivity_analysis() to see the effect of ",
            "these choices.", call. = FALSE)
  }
}

#' @export
print.model_comparison <- function(x, ...) {
  # two narrow tables, so that they do not wrap; numbers right-aligned with two
  # decimals, text left-aligned (padded, header included)
  show <- function(columns) {
    t <- x$table[columns]
    num <- vapply(t, is.double, logical(1))
    t[num] <- lapply(t[num], function(v) format(round(v, 2), nsmall = 2))
    # exact leave-one-out has no Pareto k
    if ("high_k" %in% names(t)) t$high_k <- ifelse(is.na(t$high_k), "-", t$high_k)
    text <- vapply(t, is.character, logical(1)) & !num & names(t) != "high_k"
    for (j in which(text)) {
      width <- max(nchar(c(t[[j]], names(t)[j])))
      t[[j]] <- formatC(t[[j]], width = -width)
      names(t)[j] <- formatC(names(t)[j], width = -width)
    }
    print(t, row.names = FALSE)
  }
  cat("Expected log predictive density of a new subject (leave-one-out)\n\n")
  show(c("model", "elpd_loo", "se_loo", "elpd_diff", "se_diff", "p_loo"))
  cat("\nReliability of the estimates and fit to the data\n\n")
  show(c("model", "loo_method", "high_k", "elpd_waic", "p_waic", "ppc_pvalue"))
  cat("\nelpd_diff: difference from the first model, with standard error se_diff.\n",
      "high_k: cells with Pareto k > 0.7 (- for exact leave-one-out, which has no k).\n",
      "ppc_pvalue: posterior predictive p-value of the chi-square discrepancy\n",
      "between observed and expected counts.\n", sep = "")
  invisible(x)
}

#' Plot the posterior predictive checks of a model comparison
#'
#' For each model, the observed count of each cell of the table (verified
#' diseased, verified healthy and unverified subjects of each response profile
#' of the tests) is shown with the interval of its replicated counts under the
#' posterior predictive distribution. Observed counts outside their intervals
#' are cells that the model does not reproduce. The counts are on a square root
#' scale, so that small and large cells are both visible. The label of each model
#' reports the posterior predictive p-value of the chi-square discrepancy of
#' the whole table.
#'
#' @param x Result of [compare_models()].
#' @param ... Not used.
#' @return A `ggplot` object.
#' @export
plot.model_comparison <- function(x, ...) {
  d <- x$checks
  p <- stats::setNames(x$table$ppc_pvalue, x$table$model)
  label <- function(m) sprintf("%s\n(p = %.2f)", m, p[m])
  d$model <- factor(label(d$model), levels = label(x$table$model))
  d$profile <- factor(d$profile, levels = rev(unique(d$profile)))
  ggplot2::ggplot(d, ggplot2::aes(y = .data$profile)) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = .data$lower, xmax = .data$upper,
                                         colour = "replicated"), linewidth = 1.5) +
    ggplot2::geom_point(ggplot2::aes(x = .data$observed, colour = "observed"), size = 2) +
    ggplot2::scale_colour_manual(values = c(replicated = colours[["posterior"]],
                                            observed = colours[["prior"]]), name = NULL) +
    ggplot2::scale_x_sqrt(guide = ggplot2::guide_axis(check.overlap = TRUE)) +
    ggplot2::facet_grid(rows = ggplot2::vars(.data$model), cols = ggplot2::vars(.data$cell),
                        scales = "free_x") +
    ggplot2::labs(x = "count (square root scale)", y = "profile",
                  title = "Posterior predictive check of the observed cells") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "top")
}
