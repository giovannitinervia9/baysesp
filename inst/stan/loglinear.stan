// baysesp, road 2: log-linear models of order m (specification, sections 5.2,
// 5.3 and 11.4). The distribution of the profiles in class d is
// theta_d = softmax(X_m * eta_d), where the columns of the design matrix X_m are
// the subsets of at most m tests. Order 1 is conditional independence.
//
// Verification: under MAR it is ignorable and omitted. Under MNAR,
// logit(lambda_{t,d}) = alpha_t + beta d, with one beta or one per profile,
// or with beta fixed (for the sensitivity analysis). Never-verified profiles
// have lambda = 0. Scalar parameters are on the unconstrained scale, with
// priors added by prior_lpdf(), see include/priors.stan.

functions {
#include include/priors.stan
#include include/likelihood.stan
}

data {
  int<lower=2> J;                                // number of profiles, 2^K
  int<lower=1> P;                                // parameters per class
  matrix[J, P] X;                                // design matrix X_m
  array[P] int<lower=1> order;                   // order of each column (1 = main effect)
  array[J] int<lower=0> n1;                      // verified diseased
  array[J] int<lower=0> n0;                      // verified healthy
  array[J] int<lower=0> u;                       // unverified

  int prevalence_family;                         // prior on the prevalence
  vector[3] prevalence_params;
  int main_family;                               // prior on the main effects
  vector[3] main_params;
  int interaction_family;                        // prior on the interactions
  vector[3] interaction_params;

  int<lower=0, upper=1> mnar;                    // 0 = MAR, 1 = MNAR
  int<lower=0> n_beta;                           // 0 (MAR or fixed beta), 1, or one per verified profile
  real beta_value;                               // fixed beta, used if mnar and n_beta == 0
  int alpha_family;                         // prior on alpha (used only if mnar)
  vector[3] alpha_params;
  int beta_family;                               // prior on beta (used only if n_beta > 0)
  vector[3] beta_params;
}

transformed data {
  array[J] int verified_pos = verified_positions(n1, n0);
  int J1 = max(verified_pos);                    // profiles with verified subjects
  if (n_beta > 1 && n_beta != J1) reject("n_beta must be 0, 1 or the number of verified profiles");
}

parameters {
  real logit_prevalence;
  vector[P] eta0;                                // healthy
  vector[P] eta1;                                // diseased
  vector[mnar ? J1 : 0] alpha;                   // log odds of verification of the healthy
  vector[n_beta] beta;
}

model {
  logit_prevalence ~ prior(prevalence_family, prevalence_params);
  for (p in 1:P) {
    if (order[p] == 1) {
      eta0[p] ~ prior(main_family, main_params);
      eta1[p] ~ prior(main_family, main_params);
    } else {
      eta0[p] ~ prior(interaction_family, interaction_params);
      eta1[p] ~ prior(interaction_family, interaction_params);
    }
  }
  for (a in 1:num_elements(alpha)) alpha[a] ~ prior(alpha_family, alpha_params);
  for (i in 1:n_beta) beta[i] ~ prior(beta_family, beta_params);

  target += cells_lp(logit_prevalence, log_softmax(X * eta0), log_softmax(X * eta1),
                     n1, n0, u, verified_pos, mnar, n_beta, beta_value, alpha, beta);
}

generated quantities {
  real prevalence = inv_logit(logit_prevalence);
  vector[J] theta_healthy = softmax(X * eta0);
  vector[J] theta_diseased = softmax(X * eta1);
}
