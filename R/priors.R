# Prior distributions ---------------------------------------------------------
#
# Every prior is an S7 object of a subclass of `prior_distribution`. Two
# read-only properties are shared by all of them:
# - `family`:  the distribution family ("beta", "normal", "student", "dirichlet");
# - `support`: the scale the prior acts on. "probability" priors act on a
#   probability; "real" priors act on an unconstrained parameter and, when they
#   are given to a probability, on its logit; "simplex" priors act on a vector
#   of probabilities that sum to one.

#' Prior distributions
#'
#' `prior_distribution` is the abstract parent class of all prior distributions.
#' Use the constructors [prior_beta()], [elicit_beta()], [prior_normal()],
#' [prior_student()] and [prior_dirichlet()] to create priors.
#'
#' @param family Name of the distribution family.
#' @param support Scale the prior acts on: `"probability"`, `"real"` or
#'   `"simplex"`.
#' @export
prior_distribution <- S7::new_class(
  "prior_distribution",
  abstract = TRUE,
  properties = list(
    family = S7::class_character,
    support = S7::class_character
  )
)

# helpers for validators ------------------------------------------------------

is_scalar_number <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x)
}

# beta -------------------------------------------------------------------------

#' Beta prior with mean and precision
#'
#' A Beta distribution for a probability, written in terms of its mean and its
#' precision (the sum of the two shape parameters). The usual shape parameters
#' are `shape1 = mean * precision` and `shape2 = (1 - mean) * precision`; they
#' are available as the read-only properties `@shape1` and `@shape2`.
#'
#' @param mean Mean of the distribution, in (0, 1).
#' @param precision Precision, a positive number: the larger it is, the more
#'   concentrated the distribution around its mean. `prior_beta(0.5, 2)` is the
#'   uniform distribution.
#' @return An object of class `prior_beta`.
#' @seealso [elicit_beta()] to obtain the precision from a plausible upper value.
#' @examples
#' prior_beta(mean = 0.5, precision = 2)   # uniform
#' p <- prior_beta(mean = 0.1, precision = 50)
#' p@shape1
#' p@shape2
#' @export
prior_beta <- S7::new_class(
  "prior_beta",
  parent = prior_distribution,
  properties = list(
    mean = S7::class_double,
    precision = S7::class_double,
    shape1 = S7::new_property(
      S7::class_double,
      getter = function(self) self@mean * self@precision
    ),
    shape2 = S7::new_property(
      S7::class_double,
      getter = function(self) (1 - self@mean) * self@precision
    )
  ),
  constructor = function(mean = 0.5, precision = 2) {
    S7::new_object(
      S7::S7_object(),
      family = "beta", support = "probability",
      mean = as.double(mean), precision = as.double(precision)
    )
  },
  validator = function(self) {
    if (!is_scalar_number(self@mean) || self@mean <= 0 || self@mean >= 1) {
      "@mean must be a single number in (0, 1)"
    } else if (!is_scalar_number(self@precision) || self@precision <= 0) {
      "@precision must be a single positive number"
    }
  }
)

#' Beta prior elicited from a mean and a plausible upper value
#'
#' Builds a [prior_beta()] from two simple statements about a probability:
#' its most plausible value (`mean`) and a value that it is unlikely to exceed
#' (`upper`), in the sense that `P(X <= upper) = prob`. For example, "the
#' proportion of diseased among never-verified subjects is about 0.5% and
#' almost certainly below 1.5%" is `elicit_beta(mean = 0.005, upper = 0.015)`.
#'
#' The probability `P(X <= upper)` of a Beta distribution with a given mean is
#' not monotone in the precision: there are usually two solutions, a U-shaped
#' Beta with its mass at the extremes and a regular one. The function returns
#' the regular one, the root to the right of the minimum.
#'
#' @param mean Most plausible value, in (0, 1).
#' @param upper Plausible upper value, greater than `mean` and lower than 1.
#' @param prob Probability that the value is below `upper`.
#' @return An object of class [prior_beta()].
#' @examples
#' p <- elicit_beta(mean = 0.005, upper = 0.015)
#' p
#' stats::pbeta(0.015, p@shape1, p@shape2)
#' @export
elicit_beta <- function(mean, upper, prob = 0.95) {
  if (!is_scalar_number(mean) || mean <= 0 || mean >= 1) {
    stop("`mean` must be a single number in (0, 1).", call. = FALSE)
  }
  if (!is_scalar_number(upper) || upper <= mean || upper >= 1) {
    stop("`upper` must be a single number greater than `mean` and lower than 1.",
         call. = FALSE)
  }
  if (!is_scalar_number(prob) || prob <= 0.5 || prob >= 1) {
    stop("`prob` must be a single number in (0.5, 1).", call. = FALSE)
  }
  gap <- function(log_precision) {
    precision <- exp(log_precision)
    stats::pbeta(upper, mean * precision, (1 - mean) * precision) - prob
  }
  range <- c(log(1e-3), log(1e8))
  log_min <- stats::optimize(gap, range)$minimum
  if (gap(log_min) > 0) {
    stop("No Beta distribution with mean ", mean, " has P(X <= ", upper, ") = ",
         prob, ": increase `upper` or `prob`.", call. = FALSE)
  }
  root <- stats::uniroot(gap, c(log_min, range[2]), tol = 1e-10)$root
  prior_beta(mean = mean, precision = exp(root))
}

# normal -----------------------------------------------------------------------

#' Normal prior
#'
#' A normal distribution for an unconstrained parameter. When it is given to a
#' probability, it acts on its logit.
#'
#' @param mean Mean.
#' @param sd Standard deviation, a positive number.
#' @return An object of class `prior_normal`.
#' @examples
#' prior_normal(0, 1.5)
#' @export
prior_normal <- S7::new_class(
  "prior_normal",
  parent = prior_distribution,
  properties = list(
    mean = S7::class_double,
    sd = S7::class_double
  ),
  constructor = function(mean = 0, sd = 1) {
    S7::new_object(
      S7::S7_object(),
      family = "normal", support = "real",
      mean = as.double(mean), sd = as.double(sd)
    )
  },
  validator = function(self) {
    if (!is_scalar_number(self@mean)) {
      "@mean must be a single finite number"
    } else if (!is_scalar_number(self@sd) || self@sd <= 0) {
      "@sd must be a single positive number"
    }
  }
)

# Student t --------------------------------------------------------------------

#' Student t prior
#'
#' A Student t distribution for an unconstrained parameter, with heavier tails
#' than the normal. When it is given to a probability, it acts on its logit.
#'
#' @param df Degrees of freedom, a positive number.
#' @param location Location (the centre of the distribution).
#' @param scale Scale, a positive number.
#' @return An object of class `prior_student`.
#' @examples
#' prior_student(df = 3, location = 0, scale = 1)
#' @export
prior_student <- S7::new_class(
  "prior_student",
  parent = prior_distribution,
  properties = list(
    df = S7::class_double,
    location = S7::class_double,
    scale = S7::class_double
  ),
  constructor = function(df = 3, location = 0, scale = 1) {
    S7::new_object(
      S7::S7_object(),
      family = "student", support = "real",
      df = as.double(df), location = as.double(location), scale = as.double(scale)
    )
  },
  validator = function(self) {
    if (!is_scalar_number(self@df) || self@df <= 0) {
      "@df must be a single positive number"
    } else if (!is_scalar_number(self@location)) {
      "@location must be a single finite number"
    } else if (!is_scalar_number(self@scale) || self@scale <= 0) {
      "@scale must be a single positive number"
    }
  }
)

# Dirichlet --------------------------------------------------------------------

#' Dirichlet prior
#'
#' A Dirichlet distribution for a vector of probabilities that sum to one, such
#' as the frequencies of the test result profiles.
#'
#' @param concentration Concentration parameters: a single positive number,
#'   repeated for every component, or a vector with one positive value per
#'   component. `prior_dirichlet(1)` is the uniform distribution on the simplex.
#' @return An object of class `prior_dirichlet`.
#' @examples
#' prior_dirichlet(1)
#' prior_dirichlet(c(2, 1, 1, 1))
#' @export
prior_dirichlet <- S7::new_class(
  "prior_dirichlet",
  parent = prior_distribution,
  properties = list(
    concentration = S7::class_double
  ),
  constructor = function(concentration = 1) {
    S7::new_object(
      S7::S7_object(),
      family = "dirichlet", support = "simplex",
      concentration = as.double(concentration)
    )
  },
  validator = function(self) {
    x <- self@concentration
    if (length(x) < 1L || anyNA(x) || any(!is.finite(x)) || any(x <= 0)) {
      "@concentration must contain positive finite numbers"
    }
  }
)

# printing ---------------------------------------------------------------------

format_number <- function(x) format(x, digits = 4)

prior_label <- S7::new_generic("prior_label", "prior")

S7::method(prior_label, prior_beta) <- function(prior) {
  sprintf("Beta(mean = %s, precision = %s), i.e. Beta(%s, %s)",
          format_number(prior@mean), format_number(prior@precision),
          format_number(prior@shape1), format_number(prior@shape2))
}
S7::method(prior_label, prior_normal) <- function(prior) {
  sprintf("Normal(mean = %s, sd = %s)", format_number(prior@mean), format_number(prior@sd))
}
S7::method(prior_label, prior_student) <- function(prior) {
  sprintf("Student-t(df = %s, location = %s, scale = %s)", format_number(prior@df),
          format_number(prior@location), format_number(prior@scale))
}
S7::method(prior_label, prior_dirichlet) <- function(prior) {
  sprintf("Dirichlet(concentration = %s)",
          paste(format_number(prior@concentration), collapse = ", "))
}

S7::method(print, prior_distribution) <- function(x, ...) {
  scale <- switch(x@support,
                  probability = "on a probability",
                  real = "on an unconstrained parameter (on the logit of a probability)",
                  simplex = "on a vector of probabilities")
  cat("<", sub("^baysesp::", "", class(x)[1]), "> ", prior_label(x), "\n", "  acts ", scale, "\n", sep = "")
  invisible(x)
}

# coding for Stan --------------------------------------------------------------
#
# In the Stan programs every scalar parameter is declared on the unconstrained
# scale and its prior is added by prior_lpdf() (inst/stan/include/priors.stan),
# which receives a family code and a parameter vector of fixed length 3.

stan_prior_families <- c(normal = 1L, student = 2L, beta = 3L)

prior_code <- S7::new_generic("prior_code", "prior")

S7::method(prior_code, prior_normal) <- function(prior) {
  list(family = stan_prior_families[["normal"]], params = c(prior@mean, prior@sd, 0))
}
S7::method(prior_code, prior_student) <- function(prior) {
  list(family = stan_prior_families[["student"]],
       params = c(prior@df, prior@location, prior@scale))
}
S7::method(prior_code, prior_beta) <- function(prior) {
  list(family = stan_prior_families[["beta"]], params = c(prior@shape1, prior@shape2, 0))
}
S7::method(prior_code, prior_dirichlet) <- function(prior) {
  stop("A Dirichlet prior is passed to Stan as its concentration vector, ",
       "not through prior_code().", call. = FALSE)
}

# compatibility with the parameter a prior is given to -------------------------
#
# target: the natural scale of the parameter, "probability", "real" or "simplex".
# A "real" prior can be given to a probability (it then acts on the logit); a
# "probability" prior only to a probability; a "simplex" prior only to a simplex.

check_prior <- function(prior, target, argument) {
  if (!S7::S7_inherits(prior, prior_distribution)) {
    stop("`", argument, "` must be a prior distribution, such as prior_beta() ",
         "or prior_normal().", call. = FALSE)
  }
  allowed <- switch(target,
                    probability = c("probability", "real"),
                    real = "real",
                    simplex = "simplex")
  if (!prior@support %in% allowed) {
    stop("`", argument, "` is a ", prior@family, " prior, which cannot be given ",
         "to a parameter on the ", target, " scale.", call. = FALSE)
  }
  invisible(prior)
}
