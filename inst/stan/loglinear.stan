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
  int intercepts_family;                         // prior on alpha (used only if mnar)
  vector[3] intercepts_params;
  int beta_family;                               // prior on beta (used only if n_beta > 0)
  vector[3] beta_params;
}

transformed data {
  int J1 = 0;                                    // profiles with verified subjects
  array[J] int verified_pos;                     // position among them, 0 if never verified
  for (j in 1:J) {
    if (n1[j] + n0[j] > 0) {
      J1 += 1;
      verified_pos[j] = J1;
    } else {
      verified_pos[j] = 0;
    }
  }
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
  vector[J] log_theta0 = log_softmax(X * eta0);
  vector[J] log_theta1 = log_softmax(X * eta1);
  real log_p1 = log_inv_logit(logit_prevalence);
  real log_p0 = log1m_inv_logit(logit_prevalence);

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
  for (a in 1:num_elements(alpha)) alpha[a] ~ prior(intercepts_family, intercepts_params);
  for (i in 1:n_beta) beta[i] ~ prior(beta_family, beta_params);

  // aggregated likelihood, one profile at a time (specification, section 3.4)
  for (j in 1:J) {
    real lp1 = log_p1 + log_theta1[j];           // log P(D = 1, T = t)
    real lp0 = log_p0 + log_theta0[j];           // log P(D = 0, T = t)
    int a = verified_pos[j];
    if (a == 0) {
      target += u[j] * log_sum_exp(lp1, lp0);
    } else if (mnar == 0) {
      target += n1[j] * lp1 + n0[j] * lp0 + u[j] * log_sum_exp(lp1, lp0);
    } else {
      real b = n_beta == 0 ? beta_value : beta[n_beta == 1 ? 1 : a];
      real v1 = alpha[a] + b;                    // logit lambda_{t,1}
      real v0 = alpha[a];                        // logit lambda_{t,0}
      target += n1[j] * (lp1 + log_inv_logit(v1)) + n0[j] * (lp0 + log_inv_logit(v0));
      target += u[j] * log_sum_exp(lp1 + log1m_inv_logit(v1), lp0 + log1m_inv_logit(v0));
    }
  }
}

generated quantities {
  real prevalence = inv_logit(logit_prevalence);
  vector[J] theta_healthy = softmax(X * eta0);
  vector[J] theta_diseased = softmax(X * eta1);
}
