#' @include sensitivity.R generics.R
NULL

# Prior-posterior comparison and plots -------------------------------------------
#
# All plots use {ggplot2} with theme_bw() and return the ggplot object, which
# can be modified further. Colours: posterior blue and prior orange, a pair
# checked for colour vision deficiencies, also told apart by filled and hollow
# points.

colours <- c(posterior = "#2a78d6", prior = "#eb6834")

#' Compare prior and posterior
#'
#' Fits the same model again sampling from the prior alone, and compares the
#' prior and posterior distributions of the accuracy of the tests (or of
#' rules) and, for the saturated model, of the parameters the data cannot
#' estimate: the proportion of diseased in never-verified profiles (`q`) and
#' the log odds ratio of verification (`beta`) when `verification = mnar()` (Kosinski and Barnhart, 2003).
#' This shows how much each result depends on the data and how much on the
#' assumptions: a posterior equal to its prior is not learned from the data
#' (Buzoianu and Kadane, 2008; Gustafson, 2015). For the structured models
#' it is also a prior simulation, which shows the prior that their priors induce on the
#' accuracy.
#'
#' @template accuracy
#' @inheritParams sensitivity_analysis
#' @param ... Arguments passed to [baysesp()] for the prior fit, such as
#'   `chains` and `seed`.
#' @return A data frame of class `prior_posterior` with columns `quantity`,
#'   `distribution` (`"prior"` or `"posterior"`), `mean`, `median`, `lower`
#'   and `upper`. It has a `plot()` method.
#' @references
#' Buzoianu, M. and Kadane, J. B. (2008). Adjusting for verification bias in diagnostic test evaluation: A Bayesian approach. *Statistics in Medicine*, 27(13), 2453–2473. \doi{10.1002/sim.3099}
#'
#' Gustafson, P. (2015). *Bayesian Inference for Partially Identified Models: Exploring the Limits of Limited Data*. Chapman & Hall/CRC, Boca Raton, FL.
#'
#' Kosinski, A. S. and Barnhart, H. X. (2003). Accounting for nonignorable verification bias in assessment of diagnostic tests. *Biometrics*, 59(1), 163–171. \doi{10.1111/1541-0420.00019}
#'
#' Torrance-Rynard, V. L. and Walter, S. D. (1997). Effects of dependent errors in the assessment of diagnostic test performance. *Statistics in Medicine*, 16(19), 2157–2175. \doi{10.1002/(SICI)1097-0258(19971015)16:19<2157::AID-SIM653>3.0.CO;2-X}
#'
#' Vacek, P. M. (1985). The effect of conditional dependence on the evaluation of diagnostic tests. *Biometrics*, 41(4), 959–968. \doi{10.2307/2530967}
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- simulate_study(n = 2000, prevalence = 0.2,
#'                     sensitivity = c(A = 0.9, B = 0.8),
#'                     specificity = c(A = 0.85, B = 0.95),
#'                     verification_prob = c(0, 0.9, 0.9, 1))
#' d <- baysesp_data(d, tests = c("A", "B"), verified = "verified", disease = "disease")
#' fit <- baysesp(d, saturated(unverified_prevalence_prior = elicit_beta(0.03, 0.06)),
#'                seed = 1)
#' pp <- prior_posterior(fit, seed = 1)
#' pp
#' plot(pp)
#' }
#' @export
prior_posterior <- function(fit, rules = NULL, level = 0.95, ...) {
  if (!S7::S7_inherits(fit, baysesp_fit)) {
    stop("`fit` must be a model fitted with baysesp().", call. = FALSE)
  }
  check_level(level)
  sdata <- stan_data(fit@model, fit@model@verification, fit@data)
  sdata$prior_only <- 1L
  args <- utils::modifyList(list(chains = 4, iter_warmup = 1000, iter_sampling = 1000,
                                 seed = NULL), list(...))
  prior_fit <- do.call(fit_stan, c(list(fit@model, fit@data, sdata, fit@backend), args))

  summarise <- function(f) {
    acc <- if (is.null(rules)) test_accuracy(f, level) else rule_accuracy(f, rules, level)
    acc$identified <- NULL
    m <- unclass(posterior::as_draws_matrix(f@draws))
    extra <- grep("^(q|beta)\\[", colnames(m), value = TRUE)
    if (length(extra) > 0L) {
      acc <- rbind(acc, summarise_draws(m[, extra, drop = FALSE], level, identified = NA)[, 1:5])
    }
    acc
  }
  prior <- summarise(prior_fit)
  post <- summarise(fit)
  out <- rbind(cbind(post[, 1, drop = FALSE], distribution = "prior", prior[, -1]),
               cbind(post[, 1, drop = FALSE], distribution = "posterior", post[, -1]))
  out <- out[order(match(out$quantity, post$quantity), out$distribution != "prior"), ]
  rownames(out) <- NULL
  class(out) <- c("prior_posterior", "data.frame")
  out
}

#' @export
print.prior_posterior <- function(x, ...) {
  print(as.data.frame(unclass(x)), row.names = FALSE, digits = 3, ...)
  invisible(x)
}

#' Plot the comparison of prior and posterior
#'
#' @param x Result of [prior_posterior()].
#' @param ... Not used.
#' @return A `ggplot` object.
#' @export
plot.prior_posterior <- function(x, ...) {
  x <- as.data.frame(unclass(x))
  interval_plot(x, group = "distribution")
}

# Plot of a fitted model: the result of test_accuracy().
S7::method(plot_generic, baysesp_fit) <- function(x, ..., level = 0.95) {
  interval_plot(test_accuracy(x, level)) +
    ggplot2::labs(title = paste0("Accuracy of the tests (", 100 * level, "% intervals)"))
}

# Plot of a prior: its density on the natural scale; a prior on the real line
# is also shown as the prior it induces on a probability (through the logit).
S7::method(plot_generic, prior_distribution) <- function(x, ...) {
  if (x@family == "dirichlet") {
    stop("Plot the marginals of a Dirichlet prior: each component is Beta(a_j, sum(a) - a_j).",
         call. = FALSE)
  }
  p <- seq(0.001, 0.999, length.out = 500)
  if (x@support == "probability") {
    d <- data.frame(value = p, density = stats::dbeta(p, x@shape1, x@shape2),
                    scale = "probability")
  } else {
    dens <- function(v) {
      switch(x@family,
             normal = stats::dnorm(v, x@mean, x@sd),
             student = stats::dt((v - x@location) / x@scale, x@df) / x@scale)
    }
    centre <- if (x@family == "normal") x@mean else x@location
    spread <- if (x@family == "normal") x@sd else x@scale * 3
    v <- seq(centre - 4 * spread, centre + 4 * spread, length.out = 500)
    d <- rbind(data.frame(value = v, density = dens(v), scale = "logit scale"),
               data.frame(value = p, density = dens(stats::qlogis(p)) / (p * (1 - p)),
                          scale = "induced on a probability"))
    d$scale <- factor(d$scale, levels = c("logit scale", "induced on a probability"))
  }
  ggplot2::ggplot(d, ggplot2::aes(x = .data$value, y = .data$density)) +
    ggplot2::geom_line(colour = colours[["posterior"]], linewidth = 0.8) +
    ggplot2::facet_wrap(ggplot2::vars(.data$scale), scales = "free") +
    ggplot2::labs(title = prior_label(x), x = NULL, y = "density") +
    ggplot2::theme_bw()
}

# Points with intervals, one row per quantity; with `group`, one coloured
# point and interval per group, side by side.
interval_plot <- function(x, group = NULL) {
  levels <- unique(x$quantity)
  x$y <- match(x$quantity, rev(levels))
  if (!is.null(group)) x$y <- x$y + ifelse(x[[group]] == "prior", 0.15, -0.15)
  p <- ggplot2::ggplot(x, ggplot2::aes(x = .data$median, xmin = .data$lower,
                                       xmax = .data$upper, y = .data$y))
  if (is.null(group)) {
    p <- p + ggplot2::geom_pointrange(colour = colours[["posterior"]])
  } else {
    p <- p +
      ggplot2::geom_pointrange(ggplot2::aes(colour = .data[[group]], shape = .data[[group]]),
                               fill = "white") +
      ggplot2::scale_colour_manual(values = colours, breaks = c("prior", "posterior"),
                                   name = NULL) +
      ggplot2::scale_shape_manual(values = c(prior = 21, posterior = 19),
                                  breaks = c("prior", "posterior"), name = NULL)
  }
  if (min(x$lower, na.rm = TRUE) >= 0 && max(x$upper, na.rm = TRUE) <= 1) {
    p <- p + ggplot2::scale_x_continuous(limits = c(0, 1))   # probabilities: whole scale
  }
  p +
    ggplot2::scale_y_continuous(breaks = seq_along(levels), labels = rev(levels),
                                minor_breaks = NULL) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "top")
}
