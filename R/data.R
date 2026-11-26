#' @include generics.R
NULL

# Study data --------------------------------------------------------------------
#
# The model only needs the counts of the 3J cells of the aggregated table
# (specification, section 2.3): for every response profile t, the verified
# diseased n_{t,1}, the verified healthy n_{t,0} and the unverified u_t.
# Profiles are in binary order: profile t is in row j(t) = 1 + sum_k t_k 2^(k-1),
# so the first row is "all tests negative" and the first test varies fastest.

#' Study data
#'
#' Summarises a study in which `K` binary tests are performed on every subject
#' and the disease status is ascertained by the reference standard only for the
#' verified subjects. The data are reduced, without loss of information, to the
#' counts of verified diseased, verified healthy and unverified subjects for
#' each of the `J = 2^K` response profiles.
#'
#' @param data A data frame. By default one row per subject; with `counts`, one
#'   row per combination of test results, verification and disease status.
#' @param tests Names of the columns with the test results.
#' @param verified Name of the column that says whether the subject was verified.
#' @param disease Name of the column with the disease status. It may be missing
#'   (`NA`) only for unverified subjects, and it is ignored for them.
#' @param counts Optional name of a column with the number of subjects in each
#'   row, for data that are already aggregated.
#'
#' @details Test results, verification and disease status must be coded as
#'   `0`/`1` or as logical values, with `1` (`TRUE`) meaning positive, verified
#'   and diseased.
#'
#'   The object has the properties `@tests`, `@profiles` (a `J x K` matrix of the
#'   response profiles, in binary order: the first row has all tests negative
#'   and the first test varies fastest), `@n_diseased`, `@n_healthy` and
#'   `@n_unverified` (counts for each profile) and the read-only `@never_verified`
#'   (profiles without any verified subject). `as.data.frame()` returns the
#'   aggregated table.
#'
#' @return An object of class `baysesp_data`.
#' @examples
#' d <- simulate_study(n = 1000, prevalence = 0.1,
#'                     sensitivity = c(A = 0.9, B = 0.8),
#'                     specificity = c(A = 0.95, B = 0.9),
#'                     verification_prob = c(0.1, 0.8, 0.8, 1))
#' baysesp_data(d, tests = c("A", "B"), verified = "verified", disease = "disease")
#'
#' # aggregated data
#' agg <- data.frame(T1 = c(0, 1, 0, 1, 0, 1),
#'                   V  = c(1, 1, 1, 1, 0, 0),
#'                   D  = c(0, 0, 1, 1, NA, NA),
#'                   n  = c(40, 5, 2, 30, 300, 20))
#' baysesp_data(agg, tests = "T1", verified = "V", disease = "D", counts = "n")
#' @export
baysesp_data <- S7::new_class(
  "baysesp_data",
  properties = list(
    tests = S7::class_character,
    profiles = S7::class_integer,
    n_diseased = S7::class_integer,
    n_healthy = S7::class_integer,
    n_unverified = S7::class_integer,
    never_verified = S7::new_property(
      S7::class_logical,
      getter = function(self) self@n_diseased + self@n_healthy == 0L
    )
  ),
  constructor = function(data, tests, verified, disease, counts = NULL) {
    if (!is.data.frame(data)) {
      stop("`data` must be a data frame.", call. = FALSE)
    }
    if (!is.character(tests) || length(tests) == 0L || anyDuplicated(tests)) {
      stop("`tests` must be a character vector of distinct column names.",
           call. = FALSE)
    }
    columns <- c(tests, verified, disease, counts)
    missing_columns <- setdiff(columns, names(data))
    if (length(missing_columns) > 0L) {
      stop("Columns not found in `data`: ",
           paste0("`", missing_columns, "`", collapse = ", "), ".", call. = FALSE)
    }

    test_values <- vapply(tests, function(k) as_binary(data[[k]], k), integer(nrow(data)))
    test_values <- matrix(test_values, nrow = nrow(data))
    if (anyNA(test_values)) {
      stop("Test results must not be missing.", call. = FALSE)
    }
    v <- as_binary(data[[verified]], verified)
    if (anyNA(v)) {
      stop("`", verified, "` must not be missing.", call. = FALSE)
    }
    d <- as_binary(data[[disease]], disease)
    if (anyNA(d[v == 1L])) {
      stop("`", disease, "` must not be missing for verified subjects.", call. = FALSE)
    }

    if (is.null(counts)) {
      w <- rep(1, nrow(data))
    } else {
      w <- data[[counts]]
      if (!is.numeric(w) || anyNA(w) || any(w < 0) || any(w != round(w))) {
        stop("`", counts, "` must contain non-negative whole numbers.", call. = FALSE)
      }
    }

    K <- length(tests)
    J <- 2L^K
    profile <- 1L + drop(test_values %*% 2L^(seq_len(K) - 1L))
    # cell: 1 = verified diseased, 2 = verified healthy, 3 = unverified
    cell <- ifelse(v == 0L, 3L, ifelse(d == 1L, 1L, 2L))
    table <- stats::xtabs(w ~ factor(profile, levels = seq_len(J)) +
                            factor(cell, levels = 1:3))

    profiles <- as.matrix(expand.grid(rep(list(0:1), K)))
    dimnames(profiles) <- list(NULL, tests)

    S7::new_object(
      S7::S7_object(),
      tests = tests,
      profiles = profiles,
      n_diseased = as.integer(table[, 1]),
      n_healthy = as.integer(table[, 2]),
      n_unverified = as.integer(table[, 3])
    )
  },
  validator = function(self) {
    J <- as.integer(2^length(self@tests))
    n <- list(self@n_diseased, self@n_healthy, self@n_unverified)
    if (!identical(dim(self@profiles), c(J, length(self@tests)))) {
      "@profiles must be a 2^K x K matrix"
    } else if (any(lengths(n) != J)) {
      "the counts must have one value per profile"
    } else if (anyNA(unlist(n)) || any(unlist(n) < 0L)) {
      "the counts must be non-negative"
    } else if (sum(self@n_diseased + self@n_healthy) == 0L) {
      "at least one subject must be verified"
    }
  }
)

# Converts 0/1 or logical values to integer, keeping NA; anything else is an error.
as_binary <- function(x, name) {
  if (is.logical(x)) {
    return(as.integer(x))
  }
  if (!is.numeric(x) || !all(x %in% c(0, 1, NA))) {
    stop("`", name, "` must be coded as 0/1 or as logical values.", call. = FALSE)
  }
  as.integer(x)
}

S7::method(as.data.frame, baysesp_data) <- function(x, ...) {
  data.frame(x@profiles, diseased = x@n_diseased, healthy = x@n_healthy,
             unverified = x@n_unverified, check.names = FALSE)
}

# Profiles written as "(0,1)".
profile_labels <- function(data) {
  paste0("(", apply(data@profiles, 1, paste, collapse = ","), ")")
}

S7::method(print_generic, baysesp_data) <- function(x, ...) {
  K <- length(x@tests)
  n_verified <- sum(x@n_diseased + x@n_healthy)
  never <- if (!any(x@never_verified)) "none" else
    paste(profile_labels(x)[x@never_verified], collapse = " ")
  cat("<baysesp_data> ", K, if (K == 1L) " test, " else " tests, ",
      2L^K, " profiles, ",
      sum(x@n_diseased, x@n_healthy, x@n_unverified), " subjects, ",
      n_verified, " verified\n",
      "never-verified profiles: ", never, "\n\n", sep = "")
  print(as.data.frame(x), row.names = FALSE)
  invisible(x)
}
