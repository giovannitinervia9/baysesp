#' @details
#' The accuracy is defined for a *rule* \eqn{R}, which declares a subject
#' positive or negative from the response profile \eqn{t}, the vector of the
#' results of the tests: \eqn{r_R(t) = 1} if the profile is declared positive,
#' 0 otherwise. A single test is the rule "this test is positive". The
#' quantities depend on the prevalence \eqn{\pi = P(D = 1)} and on the
#' distribution of the profiles in the diseased, \eqn{\theta_{1,t} = P(T = t \mid D = 1)},
#' and in the healthy, \eqn{\theta_{0,t} = P(T = t \mid D = 0)}:
#'
#' - **sensitivity** (`Se`):
#'   \deqn{Se_R = P(R \textrm{ positive} \mid D = 1) = \sum_t r_R(t)\, \theta_{1,t},}
#'   the proportion of diseased subjects that the rule declares positive. It
#'   is also called the true positive rate (TPR);
#' - **specificity** (`Sp`):
#'   \deqn{Sp_R = P(R \textrm{ negative} \mid D = 0) = 1 - \sum_t r_R(t)\, \theta_{0,t},}
#'   the proportion of healthy subjects that the rule declares negative;
#'   \eqn{1 - Sp_R} is the false positive rate (FPR);
#' - **positive predictive value** (`PPV`):
#'   \deqn{PPV_R = P(D = 1 \mid R \textrm{ positive}) = \frac{\pi\, Se_R}{\pi\, Se_R + (1 - \pi)(1 - Sp_R)},}
#'   the probability that a subject declared positive is diseased;
#' - **negative predictive value** (`NPV`):
#'   \deqn{NPV_R = P(D = 0 \mid R \textrm{ negative}) = \frac{(1 - \pi)\, Sp_R}{(1 - \pi)\, Sp_R + \pi\, (1 - Se_R)},}
#'   the probability that a subject declared negative is healthy.
#'
#' Sensitivity and specificity describe the rule and do not depend on the
#' prevalence; the predictive values describe what a result means for a
#' subject, and change with the prevalence of the population. For combined
#' rules, such as [any_positive()], the accuracy depends on the joint
#' distribution of the tests, not only on their individual sensitivities and
#' specificities.
#'
#' Each quantity is computed from every posterior draw of \eqn{\pi},
#' \eqn{\theta_0} and \eqn{\theta_1}. The results report its posterior mean
#' and median and an equal-tailed credible interval with probability `level`.
