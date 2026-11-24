// baysesp: saturated model in the transparent parametrization.
// Identified parameters: omega (profile frequencies), v (verified fractions),
// phi1 (proportion of diseased among the verified subjects of each profile).
// Sensitivity parameters, not informed by the data and given explicit priors:
// q (proportion of diseased in never-verified profiles) and beta (MNAR).

data {
  int<lower=1> K;
  int<lower=2> J;
  matrix[K, J] A;                                // profile matrix
  int<lower=0> R;
  matrix[R, J] Rmat;                             // rule vectors r_R, one per row

  array[J] int<lower=0> n1;                      // verified diseased
  array[J] int<lower=0> n0;                      // verified healthy
  array[J] int<lower=0> u;                       // unverified
  array[J] int<lower=0, upper=1> never_verified; // 1 if the profile is never verified

  int<lower=0, upper=1> mnar;                    // 0 = MAR, 1 = MNAR with parameter beta
  vector<lower=0>[J] a_q;                        // Beta(a_q, b_q) priors on q_t
  vector<lower=0>[J] b_q;                        //   (used only for never-verified profiles)
  real beta_mean;
  real<lower=0> beta_sd;
}

transformed data {
  int J0 = sum(never_verified);                  // number of never-verified profiles
  int J1 = J - J0;
  array[J1] int idx_v;                           // indices of verified profiles
  array[J0] int idx_0;                           // indices of never-verified profiles
  array[J] int Nt;
  {
    int a = 1;
    int b = 1;
    for (j in 1:J) {
      Nt[j] = n1[j] + n0[j] + u[j];
      if (never_verified[j] == 1) { idx_0[b] = j; b += 1; }
      else { idx_v[a] = j; a += 1; }
    }
  }
}

parameters {
  simplex[J] omega;                              // P(T = t)
  vector<lower=0, upper=1>[J1] v;                // P(V = 1 | T = t), verified profiles
  vector<lower=0, upper=1>[J1] phi1;             // P(D = 1 | T = t, V = 1)
  vector<lower=0, upper=1>[J0] q;                // P(D = 1 | T = t), never-verified profiles
  vector[mnar] beta;                             // log odds ratio of verification, diseased vs healthy
}

model {
  // priors on the identified parameters: weak, the likelihood dominates
  omega ~ dirichlet(rep_vector(1, J));
  v ~ beta(1, 1);
  phi1 ~ beta(1, 1);
  // priors on the sensitivity parameters: they alone determine their posteriors
  for (b in 1:J0) q[b] ~ beta(a_q[idx_0[b]], b_q[idx_0[b]]);
  if (mnar) beta ~ normal(beta_mean, beta_sd);

  // likelihood: multinomial over profiles and independent binomials
  Nt ~ multinomial(omega);
  for (a in 1:J1) {
    int j = idx_v[a];
    n1[j] + n0[j] ~ binomial(Nt[j], v[a]);
    n1[j] ~ binomial(n1[j] + n0[j], phi1[a]);
  }
}

generated quantities {
  vector[J] psi;                                 // P(D = 1 | T = t)
  real prev;
  vector[J] theta1;
  vector[J] theta0;
  for (a in 1:J1) {
    int j = idx_v[a];
    real phi0 = mnar ? inv_logit(logit(phi1[a]) - beta[1]) : phi1[a];
    psi[j] = v[a] * phi1[a] + (1 - v[a]) * phi0;
  }
  for (b in 1:J0) psi[idx_0[b]] = q[b];
  prev = dot_product(omega, psi);
  theta1 = (omega .* psi) / prev;
  theta0 = (omega .* (1 - psi)) / (1 - prev);

  vector[K] Se = A * theta1;
  vector[K] Sp = 1 - A * theta0;
  vector[R] TPR = Rmat * theta1;
  vector[R] FPR = Rmat * theta0;
}
