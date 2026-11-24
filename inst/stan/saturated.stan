// baysesp, road 1: saturated model in the transparent parametrization
// (specification, sections 7.2 and 11.4).
//
// Identified parameters, informed by the data:
//   omega  profile frequencies, P(T = t);
//   v      verified fraction of each verified profile, P(V = 1 | T = t);
//   phi1   proportion of diseased among the verified, P(D = 1 | T = t, V = 1).
// Sensitivity parameters, absent from the likelihood (posterior = prior):
//   q      proportion of diseased in the never-verified profiles;
//   beta   log odds ratio of verification, diseased vs healthy (MNAR only).
// Scalar parameters are declared on the unconstrained (logit) scale and their
// priors are added by prior_lpdf(), see include/priors.stan.

functions {
#include include/priors.stan
}

data {
  int<lower=1> J;                                // number of profiles, 2^K
  array[J] int<lower=0> n1;                      // verified diseased
  array[J] int<lower=0> n0;                      // verified healthy
  array[J] int<lower=0> u;                       // unverified

  vector<lower=0>[J] frequency_concentration;    // Dirichlet prior on omega
  int verified_family;                           // prior on v
  vector[3] verified_params;
  int prevalence_family;                         // prior on phi1
  vector[3] prevalence_params;
  int<lower=0> J0;                               // number of never-verified profiles
  array[J0] int q_family;                        // prior on q, one per never-verified profile
  array[J0] vector[3] q_params;

  int<lower=0, upper=1> mnar;                    // 0 = MAR, 1 = MNAR with parameter beta
  int beta_family;                               // prior on beta (used only if mnar)
  vector[3] beta_params;
}

transformed data {
  int J1 = J - J0;
  array[J1] int verified_idx;                    // profiles with at least one verified subject
  array[J0] int never_idx;                       // never-verified profiles
  array[J] int N_profile;
  {
    int a = 0;
    int b = 0;
    for (j in 1:J) {
      N_profile[j] = n1[j] + n0[j] + u[j];
      if (n1[j] + n0[j] > 0) {
        a += 1;
        verified_idx[a] = j;
      } else {
        b += 1;
        never_idx[b] = j;
      }
    }
    if (b != J0) reject("J0 does not match the number of never-verified profiles");
  }
}

parameters {
  simplex[J] omega;
  vector[J1] logit_v;
  vector[J1] logit_phi1;
  vector[J0] logit_q;
  vector[mnar] beta;
}

model {
  omega ~ dirichlet(frequency_concentration);
  for (a in 1:J1) {
    logit_v[a] ~ prior(verified_family, verified_params);
    logit_phi1[a] ~ prior(prevalence_family, prevalence_params);
  }
  for (b in 1:J0) logit_q[b] ~ prior(q_family[b], q_params[b]);
  if (mnar) beta[1] ~ prior(beta_family, beta_params);

  // likelihood: multinomial over profiles and independent binomials
  N_profile ~ multinomial(omega);
  for (a in 1:J1) {
    int j = verified_idx[a];
    n1[j] + n0[j] ~ binomial_logit(N_profile[j], logit_v[a]);
    n1[j] ~ binomial_logit(n1[j] + n0[j], logit_phi1[a]);
  }
}

generated quantities {
  real prevalence;
  vector[J] theta_healthy;
  vector[J] theta_diseased;
  vector[J0] q = inv_logit(logit_q);
  {
    vector[J] psi;                               // P(D = 1 | T = t)
    for (a in 1:J1) {
      // proportion of diseased among the unverified (specification, Proposition 7.2)
      real logit_phi0 = mnar ? logit_phi1[a] - beta[1] : logit_phi1[a];
      real v = inv_logit(logit_v[a]);
      psi[verified_idx[a]] = v * inv_logit(logit_phi1[a]) + (1 - v) * inv_logit(logit_phi0);
    }
    for (b in 1:J0) psi[never_idx[b]] = q[b];
    prevalence = dot_product(omega, psi);
    theta_diseased = omega .* psi / prevalence;
    theta_healthy = omega .* (1 - psi) / (1 - prevalence);
  }
}
