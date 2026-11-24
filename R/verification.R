#' Verification mechanisms
#'
#' `verification_mechanism` is the abstract parent class of the assumptions on
#' how subjects are selected for verification by the reference standard.
#' `mar()` assumes that verification depends only on the test results
#' (missing at random): among subjects with the same profile, verified and
#' unverified subjects have the same proportion of diseased.
#'
#' @return An object of class `mar`.
#' @examples
#' mar()
#' @export
verification_mechanism <- S7::new_class("verification_mechanism", abstract = TRUE)

#' @rdname verification_mechanism
#' @export
mar <- S7::new_class("mar", parent = verification_mechanism)

S7::method(print, mar) <- function(x, ...) {
  cat("<mar> verification depends only on the test results\n")
  invisible(x)
}
