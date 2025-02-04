# data
# s1 = #(verification = 1, test = 1, disease = 1) true positive
# s2 = #(verification = 1, test = 0, disease = 1) false negative
# r1 = #(verification = 1, test = 1, disease = 0) false positive
# r2 = #(verification = 1, test = 0, disease = 0) true negative
# u1 = #(verification = 0, test = 1, disease = NA)
# u2 = #(verification = 0, test = 0, disease = NA)


# parameters
# se = sensitivity
# sp = specificity
# p = prevalence
# l11 = P(V = 1|T = 1, D = 1)
# l21 = P(V = 1|T = 1, D = 0)
# l12 = P(V = 1|T = 0, D = 1)
# l22 = P(V = 1|T = 0, D = 0)

# baysesp



#-------------------------------------------------------------------------------

#' Simulate Diagnostic Test Data with Verification Bias
#'
#' This function generates simulated diagnostic test data incorporating
#' verification bias. It uses multinomial sampling to create data based on
#' specified sensitivity, specificity, prevalence, and verification probabilities.
#'
#' @param n Integer. The total number of observations to simulate. Default is 100.
#' @param nsim Integer. The number of simulation replicates. Default is 1.
#' @param p Numeric. The prevalence of the disease in the population. Default is 0.01.
#' @param se Numeric. The sensitivity of the diagnostic test, defined as
#' \eqn{P(T = 1 | D = 1)}.
#' @param sp Numeric. The specificity of the diagnostic test, defined as
#' \eqn{P(T = 0 | D = 0)}.
#' @param l11 Numeric. Probability of verification given a positive test result
#' and a positive disease status, defined as \eqn{P(V = 1 | T = 1, D = 1)}.
#' @param l21 Numeric. Probability of verification given a positive test result
#' and a negative disease status, defined as \eqn{P(V = 1 | T = 1, D = 0)}.
#' @param l12 Numeric. Probability of verification given a negative test result
#' and a positive disease status, defined as \eqn{P(V = 1 | T = 0, D = 1)}.
#' @param l22 Numeric. Probability of verification given a negative test result
#' and a negative disease status, defined as \eqn{P(V = 1 | T = 0, D = 0)}.
#'
#' @return A matrix of counts with 6 rows and `nsim` columns.
#' The rows represent:
#' - `"s1"`: Verified true positives \eqn{(V = 1, T = 1, D = 1)}
#' - `"s2"`: Verified false negatives \eqn{(V = 1, T = 0, D = 1)}
#' - `"r1"`: Verified false positives \eqn{(V = 1, T = 1, D = 0)}
#' - `"r2"`: Verified true negatives \eqn{(V = 1, T = 0, D = 0)}
#' - `"u1"`: Unverified test positives \eqn{(V = 0, T = 1, D = NA)}
#' - `"u2"`: Unverified test negatives \eqn{(V = 0, T = 0, D = NA)}
#'
#' @examples
#' # Simulate a single dataset with default parameters
#' set.seed(123)
#' baysesp_simul()
#'
#' # Simulate 10 datasets with 500 observations each
#' set.seed(123)
#' baysesp_simul(n = 500, nsim = 10, p = 0.05, se = 0.85, sp = 0.95)
#'
#' @importFrom stats rmultinom
#' @export
baysesp_simul <- function(n = 100,
                          nsim = 1,
                          p = 0.01,
                          se = 0.8,
                          sp = 0.9,
                          l11 = 0.8,
                          l21 = 0.3,
                          l12 = 0.4,
                          l22 = 0.5) {
  p_v1_d1_t1 <- l11 * se * p                                 # s1
  p_v1_d1_t0 <- l12 * (1 - se) * p                           # s2
  p_v1_d0_t1 <- l21 * (1 - sp) * (1 - p)                     # r1
  p_v1_d0_t0 <- l22 * sp * (1 - p)                           # r2
  p_v0_t1 <- (1 - l11) * se * p + (1 - l21) * (1 - sp) * (1 - p) # u1
  p_v0_t0 <- (1 - l12) * (1 - se) * p + (1 - l22) * sp * (1 - p) # u2

  probs <- c(p_v1_d1_t1,
             p_v1_d1_t0,
             p_v1_d0_t1,
             p_v1_d0_t0,
             p_v0_t1,
             p_v0_t0)

  y <- rmultinom(n = nsim, size = n, prob = probs)

  rownames(y) <- c("s1", "s2", "r1", "r2", "u1", "u2")
  colnames(y) <- paste0("nsim", 1:nsim)
  y
}



#-------------------------------------------------------------------------------

#' Convert a Count Matrix to a Named List for `baysesp()`
#'
#' This function converts a matrix containing diagnostic test data into a named
#' list format suitable for use in `baysesp()`. It extracts the data for a specified
#' column index and returns it as a list with named elements corresponding to different
#' diagnostic test outcome categories.
#'
#' @param y Matrix. A numeric matrix where each row represents a count category,
#' and each column corresponds to a different dataset (e.g., different simulations or
#' different data sources). The matrix must have exactly six rows
#' with the following names:
#' \describe{
#'   \item{\code{"s1"}}{Count of verified true positives \eqn{(V = 1, T = 1, D = 1)}.}
#'   \item{\code{"s2"}}{Count of verified false negatives \eqn{(V = 1, T = 0, D = 1)}.}
#'   \item{\code{"r1"}}{Count of verified false positives \eqn{(V = 1, T = 1, D = 0)}.}
#'   \item{\code{"r2"}}{Count of verified true negatives \eqn{(V = 1, T = 0, D = 0)}.}
#'   \item{\code{"u1"}}{Count of unverified individuals who tested positive \eqn{(V = 0, T = 1, D = \text{NA})}.}
#'   \item{\code{"u2"}}{Count of unverified individuals who tested negative \eqn{(V = 0, T = 0, D = \text{NA})}.}
#' }
#'
#' The function assumes that row names are assigned to ensure correct mapping of
#' elements. If row names are missing, it will assign default names.
#'
#' @param index Integer. The column index in \code{y} to be extracted and converted
#' into a named list. Default is 1.
#'
#' @return A named list containing the counts for each category. If row names are
#' missing, the function assigns default names ("s1", "s2", "r1", "r2", "u1", "u2")
#' in order.
#'
#' @examples
#' # Example count matrix with named rows
#' y <- matrix(c(50, 30, 20, 40, 60, 100), ncol = 1,
#'             dimnames = list(c("s1", "s2", "r1", "r2", "u1", "u2"), "dataset1"))
#'
#' # Convert the first column to a Stan-compatible list
#' stan_data <- y_to_stan_data(y, index = 1)
#'
#' @importFrom stats setNames
#' @export
y_to_stan_data <- function(y, index = 1) {
  as.list(setNames(y[, index], rownames(y)))
}



#-------------------------------------------------------------------------------
#' Convert a Count Matrix into an Individual-Level Data Frame
#'
#' This function transforms a matrix containing counts of individuals categorized
#' by verification status, test result, and disease status into an expanded data
#' frame where each row represents an individual.
#'
#' @param y Matrix. A matrix where each row corresponds to a specific category of
#' individuals, and each column represents a different dataset (e.g., different
#' simulations or different data sources). The matrix must have exactly six rows,
#' named as follows:
#' \describe{
#'   \item{\code{"s1"}}{Count of verified true positives \eqn{(V = 1, T = 1, D = 1)}.}
#'   \item{\code{"s2"}}{Count of verified false negatives \eqn{(V = 1, T = 0, D = 1)}.}
#'   \item{\code{"r1"}}{Count of verified false positives \eqn{(V = 1, T = 1, D = 0)}.}
#'   \item{\code{"r2"}}{Count of verified true negatives \eqn{(V = 1, T = 0, D = 0)}.}
#'   \item{\code{"u1"}}{Count of unverified individuals who tested positive \eqn{(V = 0, T = 1, D = \text{NA})}.}
#'   \item{\code{"u2"}}{Count of unverified individuals who tested negative \eqn{(V = 0, T = 0, D = \text{NA})}.}
#' }
#' The column(s) of \code{y} can contain multiple datasets, where each column
#' corresponds to a different data instance (e.g., different simulations or subsets).
#'
#' @param index Integer. The column index in \code{y} to be converted into an
#' individual-level data frame. Default is 1.
#'
#' @return A data frame with three columns:
#' \describe{
#'   \item{\code{verification}}{Binary indicator: \eqn{1} if the individual was verified, \eqn{0} otherwise.}
#'   \item{\code{test}}{Binary indicator: \eqn{1} if the diagnostic test was positive, \eqn{0} if negative.}
#'   \item{\code{disease}}{Binary indicator: \eqn{1} if the individual has the disease, \eqn{0} if not, \code{NA} if unverified.}
#' }
#'
#' Each row corresponds to an individual, and the total number of rows matches the
#' sum of counts in the selected column of \code{y}.
#'
#' @examples
#' # Example with a count matrix
#' y <- matrix(c(50, 30, 20, 40, 60, 100), ncol = 1,
#'             dimnames = list(c("s1", "s2", "r1", "r2", "u1", "u2"), "dataset1"))
#'
#' # Convert the data into an individual-level format
#' df <- y_to_data_frame(y, index = 1)
#'
#' @export
y_to_data_frame <- function(y, index = 1) {
  # Define the mapping for each category
  scheme <- data.frame(
    verification = c(1, 1, 1, 1, 0, 0),
    test         = c(1, 0, 1, 0, 1, 0),
    disease      = c(1, 1, 0, 0, NA, NA),
    count        = y[, index],
    # Extract the counts from matrix y
    row.names    = c("s1", "s2", "r1", "r2", "u1", "u2")
  )

  # Expand into individual rows
  df <- scheme[rep(rownames(scheme), scheme$count), c("verification", "test", "disease")]

  # Reset row names
  rownames(df) <- NULL

  df

}



#-------------------------------------------------------------------------------

#' Convert Individual-Level Data Frame to Count-Based List for Stan
#'
#' This function processes a data frame containing individual-level diagnostic test data
#' and converts it into a named list of counts, structured for use in Bayesian modeling
#' (e.g., with Stan). The input data frame must include three binary columns:
#' \code{verification}, \code{test}, and \code{disease}.
#'
#' @param data Data frame. A dataset where each row represents an individual and includes
#' the following binary-coded columns:
#' \describe{
#'   \item{\code{verification}}{Indicator for whether the disease status was verified (\eqn{1} = verified, \eqn{0} = not verified).}
#'   \item{\code{test}}{Indicator for test result (\eqn{1} = positive, \eqn{0} = negative).}
#'   \item{\code{disease}}{Indicator for disease presence (\eqn{1} = has disease, \eqn{0} = no disease).}
#' }
#'
#' @return A named list containing counts of individuals in different verification, test, and disease status categories:
#' \describe{
#'   \item{\code{"s1"}}{Count of verified true positives \eqn{(V = 1, T = 1, D = 1)}.}
#'   \item{\code{"s2"}}{Count of verified false negatives \eqn{(V = 1, T = 0, D = 1)}.}
#'   \item{\code{"r1"}}{Count of verified false positives \eqn{(V = 1, T = 1, D = 0)}.}
#'   \item{\code{"r2"}}{Count of verified true negatives \eqn{(V = 1, T = 0, D = 0)}.}
#'   \item{\code{"u1"}}{Count of unverified individuals who tested positive \eqn{(V = 0, T = 1, D = \text{NA})}.}
#'   \item{\code{"u2"}}{Count of unverified individuals who tested negative \eqn{(V = 0, T = 0, D = \text{NA})}.}
#' }
#'
#' @examples
#' # Example dataset
#' data <- data.frame(
#'   verification = c(1, 1, 1, 1, 0, 0, 1, 0, 1, 0),
#'   test         = c(1, 0, 1, 0, 1, 0, 1, 1, 0, 0),
#'   disease      = c(1, 1, 0, 0, NA, NA, 1, NA, 0, NA)
#' )
#'
#' # Convert to Stan-compatible count list
#' stan_data <- data_frame_to_stan_data(data)
#' @importFrom dplyr filter
#' @export
data_frame_to_stan_data <- function(data) {
  test <- data$test
  verification <- data$verification
  disease = data$disease
  s1 <- data |> dplyr::filter(test == 1 &
                                verification == 1 &
                                disease == 1) |> nrow()


  s2 <- data |> dplyr::filter(test == 0 &
                                verification == 1 &
                                disease == 1) |> nrow()


  r1 <- data |> dplyr::filter(test == 1 &
                                verification == 1 &
                                disease == 0) |> nrow()


  r2 <- data |> dplyr::filter(test == 0 &
                                verification == 1 &
                                disease == 0) |> nrow()

  u1 <- data |> dplyr::filter(test == 1 &
                                verification == 0) |> nrow()

  # test negative, verification not done
  u2 <- data |> dplyr::filter(test == 0 &
                                verification == 0) |> nrow()



  list(
    s1 = s1,
    s2 = s2,
    r1 = r1,
    r2 = r2,
    u1 = u1,
    u2 = u2
  )

}


#-------------------------------------------------------------------------------

#' Generate a Stan Model Code for Bayesian Estimation
#'
#' This function generates a Stan model code for Bayesian estimation of diagnostic
#' test parameters under verification bias. The user provides a list of prior
#' distributions for each parameter, and the function returns the corresponding
#' Stan code.
#'
#' @param priors Named list. A list where each element is a character string
#' specifying the prior distribution for a parameter. The list must include
#' the following named elements:
#' \describe{
#'   \item{\code{"p"}}{Prior for disease prevalence \eqn{p}.}
#'   \item{\code{"se"}}{Prior for test sensitivity \eqn{Se = P(T = 1 | D = 1)}.}
#'   \item{\code{"sp"}}{Prior for test specificity \eqn{Sp = P(T = 0 | D = 0)}.}
#'   \item{\code{"l11"}}{Prior for verification probability given \eqn{(T = 1, D = 1)} \eqn{P(V = 1 | T = 1, D = 1)}.}
#'   \item{\code{"l21"}}{Prior for verification probability given \eqn{(T = 1, D = 0)} \eqn{P(V = 1 | T = 1, D = 0)}.}
#'   \item{\code{"l12"}}{Prior for verification probability given \eqn{(T = 0, D = 1)} \eqn{P(V = 1 | T = 0, D = 1)}.}
#'   \item{\code{"l22"}}{Prior for verification probability given \eqn{(T = 0, D = 0)} \eqn{P(V = 1 | T = 0, D = 0)}.}
#' }
#'
#' Each prior should be specified in Stan syntax (e.g., \code{"beta(2,2)"} or \code{"normal(0.5, 0.1)"}).
#'
#' @return A character string containing the complete Stan model code, including:
#' \itemize{
#'   \item \strong{Data block}: Defines the observed counts of test outcomes (\code{s1, s2, r1, r2, u1, u2}).
#'   \item \strong{Transformed Data block}: Stores the counts in a vector.
#'   \item \strong{Parameters block}: Defines the unknown parameters (\code{p, se, sp, l11, l21, l12, l22}).
#'   \item \strong{Model block}: Assigns prior distributions and specifies the likelihood using a multinomial distribution.
#' }
#'
#' @examples
#' # Define prior distributions
#' priors <- list(
#'   p = "beta(2, 2)",
#'   se = "beta(2, 2)",
#'   sp = "beta(2, 2)",
#'   l11 = "beta(2, 2)",
#'   l21 = "beta(2, 2)",
#'   l12 = "beta(2, 2)",
#'   l22 = "beta(2, 2)"
#' )
#'
#' # Generate the Stan model code
#' stan_code <- generate_stan_model(priors)
#' cat(stan_code)
#' @export
generate_stan_model <- function(priors) {
  stan_code <- paste0(
    "
  data {
  int<lower=0> s1;  // #(V = 1, T = 1, D = 1) true positive
  int<lower=0> s2;  // #(V = 1, T = 0, D = 1) false negative
  int<lower=0> r1;  // #(V = 1, T = 1, D = 0) false positive
  int<lower=0> r2;  // #(V = 1, T = 0, D = 0) true negative
  int<lower=0> u1;  // #(V = 0, T = 1)
  int<lower=0> u2;  // #(V = 0, T = 0)
}

  transformed data {
  int y[6] = {s1, s2, r1, r2, u1, u2};  // Vector of counts
}

parameters {
  real<lower=0, upper=1> p;     // Prevalence
  real<lower=0, upper=1> se;    // Sensitivity
  real<lower=0, upper=1> sp;    // Specificity
  real<lower=0, upper=1> l11;   // P(V = 1 | T = 1, D = 1)
  real<lower=0, upper=1> l21;   // P(V = 1 | T = 1, D = 0)
  real<lower=0, upper=1> l12;   // P(V = 1 | T = 0, D = 1)
  real<lower=0, upper=1> l22;   // P(V = 1 | T = 0, D = 0)
}



  model {
    p ~ ", priors[["p"]], ";
    se ~ ", priors[["se"]], ";
    sp ~ ", priors[["sp"]], ";
    l11 ~ ", priors[["l11"]], ";
    l21 ~ ", priors[["l21"]], ";
    l12 ~ ", priors[["l12"]], ";
    l22 ~ ", priors[["l22"]], ";

    vector[6] p_y;  // Probability vector for the multinomial

  p_y[1] = l11 * se * p;                          // p_s1
  p_y[2] = l12 * (1 - se) * p;                    // p_s2
  p_y[3] = l21 * (1 - sp) * (1 - p);              // p_r1
  p_y[4] = l22 * sp * (1 - p);                    // p_r2
  p_y[5] = (1 - l11) * se * p + (1 - l21) * (1 - sp) * (1 - p);  // p_u1
  p_y[6] = (1 - l12) * (1 - se) * p + (1 - l22) * sp * (1 - p);  // p_u2

  // Ensure numerical stability
  p_y = p_y / sum(p_y);


  // Multinomial likelihood
  y ~ multinomial(p_y);
  }")

  stan_code
}



#-------------------------------------------------------------------------------

#' Bayesian Estimation of Sensitivity and Specificity with Verification Bias
#'
#' This function estimates the sensitivity and specificity of a diagnostic test
#' in the presence of verification bias using a Bayesian approach implemented in Stan.
#' It constructs a Stan model based on user-specified priors, compiles and runs the model,
#' and returns the fitted Stan object.
#'
#' @param stan_data Named list. The data to be passed to Stan, typically generated
#' using \code{\link{data_frame_to_stan_data}} or \code{\link{y_to_stan_data}}.
#' It must contain the following named elements:
#' \describe{
#'   \item{\code{"s1"}}{Count of verified true positives \eqn{(V = 1, T = 1, D = 1)}.}
#'   \item{\code{"s2"}}{Count of verified false negatives \eqn{(V = 1, T = 0, D = 1)}.}
#'   \item{\code{"r1"}}{Count of verified false positives \eqn{(V = 1, T = 1, D = 0)}.}
#'   \item{\code{"r2"}}{Count of verified true negatives \eqn{(V = 1, T = 0, D = 0)}.}
#'   \item{\code{"u1"}}{Count of unverified individuals who tested positive \eqn{(V = 0, T = 1, D = \text{NA})}.}
#'   \item{\code{"u2"}}{Count of unverified individuals who tested negative \eqn{(V = 0, T = 0, D = \text{NA})}.}
#' }
#'
#' @param chains Integer. Number of Markov chains to run in the Stan model. Default is 4.
#' @param iter Integer. Total number of iterations per chain. Default is 2000.
#' @param warmup Integer. Number of warmup (burn-in) iterations per chain. Default is \code{iter/2}.
#' @param cores Integer. Number of CPU cores to use for parallel computation. Defaults to \code{getOption("mc.cores", 1L)}.
#' @param priors Named list. A list specifying the prior distributions for each parameter.
#' The list must include the following named elements:
#' \describe{
#'   \item{\code{"p"}}{Prior for disease prevalence \eqn{p}.}
#'   \item{\code{"se"}}{Prior for sensitivity \eqn{Se = P(T = 1 | D = 1)}.}
#'   \item{\code{"sp"}}{Prior for specificity \eqn{Sp = P(T = 0 | D = 0)}.}
#'   \item{\code{"l11"}}{Prior for verification probability given \eqn{(T = 1, D = 1)} \eqn{P(V = 1 | T = 1, D = 1)}.}
#'   \item{\code{"l21"}}{Prior for verification probability given \eqn{(T = 1, D = 0)} \eqn{P(V = 1 | T = 1, D = 0)}.}
#'   \item{\code{"l12"}}{Prior for verification probability given \eqn{(T = 0, D = 1)} \eqn{P(V = 1 | T = 0, D = 1)}.}
#'   \item{\code{"l22"}}{Prior for verification probability given \eqn{(T = 0, D = 0)} \eqn{P(V = 1 | T = 0, D = 0)}.}
#' }
#' Each prior should be specified using Stan's syntax (e.g., \code{"beta(1,1)"}).
#'
#' @param ... Additional arguments to pass to the \code{\link[rstan]{stan}} function.
#'
#' @return A \code{\link[rstan]{stanfit}} object containing the posterior samples
#' and summary statistics from the Bayesian estimation.
#'
#' @examples
#' # Generate example data
#' set.seed(123)
#' y <- baysesp_simul(n = 500, nsim = 1)
#' stan_data <- y_to_stan_data(y)
#'
#' # Run Bayesian estimation with default priors
#' \dontrun{
#' fit <- baysesp(stan_data)
#' print(fit)
#' }
#' @importFrom rstan stan
#' @export
baysesp <- function(stan_data,
                    chains = 4,
                    iter = 2000,
                    warmup = floor(iter / 2),
                    cores = getOption("mc.cores", 1L),
                    priors = list(
                      p = "beta(1, 1)",
                      se = "beta(1, 1)",
                      sp = "beta(1, 1)",
                      l11 = "beta(1, 1)",
                      l21 = "beta(1, 1)",
                      l12 = "beta(1, 1)",
                      l22 = "beta(1, 1)"
                    ),
                    ...) {
  model_code <- generate_stan_model(priors)

  stan(
    model_code = model_code,
    data = stan_data,
    chains = chains,
    iter = iter,
    warmup = warmup,
    cores = cores,
    ...
  )

}



#------------------------------------------------------------------------------
#' Compute Beta Distribution Parameters from Mean and Standard Deviation
#'
#' This function calculates the shape parameters (\eqn{\alpha} and \eqn{\beta})
#' of a Beta distribution given a specified mean (\eqn{\mu}) and standard deviation (\eqn{\sigma}).
#' Optionally, it also plots the resulting Beta distribution.
#'
#' @param mu Numeric. The mean of the Beta distribution. Default is 0.5.
#' @param sd Numeric. The standard deviation of the Beta distribution. Default is 0.5.
#' @param plot Logical. If \code{TRUE}, returns a plot of the Beta density function. Default is \code{TRUE}.
#'
#' @return If \code{plot = FALSE}, returns a named numeric vector with elements:
#'   \item{alpha}{Shape parameter \eqn{\alpha} of the Beta distribution.}
#'   \item{beta}{Shape parameter \eqn{\beta} of the Beta distribution.}
#'
#'   If \code{plot = TRUE}, returns a list containing:
#'   \item{par}{A named numeric vector with the computed \eqn{\alpha} and \eqn{\beta} values.}
#'   \item{plot}{A ggplot object displaying the Beta density function.}
#'
#' @examples
#' # Compute Beta parameters without plotting
#' prior_beta(mu = 0.7, sd = 0.1, plot = FALSE)
#'
#' # Compute and visualize Beta distribution
#' res <- prior_beta(mu = 0.3, sd = 0.05, plot = TRUE)
#' res$plot  # Display the plot
#'
#' @import ggplot2
#' @importFrom stats dbeta
#' @export
prior_beta <- function (mu = 0.5,
                        sd = 0.5,
                        plot = TRUE)
{
  var <- sd ^ 2
  precision <- 1 / var
  phi <- precision
  theta <- mu / (1 - mu)
  beta <- (phi * theta + (1 + theta) ^ 2) / (1 + theta) ^ 3
  alpha <- beta * theta
  par <- c(alpha = alpha, beta = beta)
  string <- paste0("beta(", par[1], ", ", par[2], ")")
  res <- list(par = par, string = string)
  if (plot) {
    x <- seq(0.001, 0.999, l = 100)
    f <- dbeta(x, alpha, beta)
    d <- data.frame(x = x, f = f)
    p <- ggplot2::ggplot(d, aes(x = x, y = f)) + ggplot2::geom_line() +
      ggplot2::labs(
        title = bquote(alpha == .(round(alpha, 3)) ~ "," ~ beta == .(round(beta, 3))),
        x = "x",
        y = expression(f(x ~ ";" ~ alpha, beta))
      ) + ggplot2::theme_bw()
    res$plot <- p
  }
  res
}
