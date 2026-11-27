// baysesp, road 2: random effects model (specification, section 5.4). Given
// the disease class d and a latent Z ~ N(0, 1), the tests are independent with
// P(T_k = 1 | D = d, Z = z) = Phi(a_dk + b_dk z). The distribution of the
// profiles integrates Z out with Gauss-Hermite quadrature:
//   theta_{d,t} = sum_q w_q prod_k Phi(a_dk + b_dk z_q)^t_k (1 - Phi(...))^(1 - t_k).
// Identifiability: changing the sign of Z and of all the loadings of a class
// leaves the model unchanged. With positive_loadings = 1 (default) all the
// loadings are positive, so only positive conditional dependence is allowed;
// with 0 only the first loading is positive, as in the specification, and the
// posterior can have a reflected mode (first loading near 0, the others with
// the opposite sign) that traps the chains. The intercepts are written through the marginal probits,
// a_dk = m_dk sqrt(1 + b_dk^2), with P(T_k = 1 | D = d) = Phi(m_dk): this
// removes most divergences when the loadings are close to 0. Verification as
// in loglinear.stan.

functions {
#include include/priors.stan
#include include/likelihood.stan

  // log theta_d for all profiles, on the log scale for numerical stability
  vector log_profiles(matrix profiles, vector a, vector b, vector nodes,
                      vector log_weights) {
    int J = rows(profiles);
    int K = cols(profiles);
    int Q = num_elements(nodes);
    matrix[Q, K] log_pos;                        // log Phi(a_k + b_k z_q)
    matrix[Q, K] log_neg;                        // log (1 - Phi(a_k + b_k z_q))
    vector[J] out;
    for (q in 1:Q) {
      for (k in 1:K) {
        // log Phi(-x) rather than log(1 - Phi(x)): finite in the far tails,
        // which the outer quadrature nodes reach
        real x = a[k] + b[k] * nodes[q];
        log_pos[q, k] = std_normal_lcdf(x);
        log_neg[q, k] = std_normal_lcdf(-x);
      }
    }
    for (j in 1:J) {
      vector[Q] terms = log_weights;
      for (k in 1:K) {
        terms += profiles[j, k] == 1 ? log_pos[, k] : log_neg[, k];
      }
      out[j] = log_sum_exp(terms);
    }
    return out;
  }
}

data {
  int<lower=0, upper=1> prior_only;             // 1 = sample from the prior only
  int<lower=2> J;                                // number of profiles, 2^K
  int<lower=2> K;                                // number of tests
  int<lower=0, upper=1> positive_loadings;       // 1 = all loadings positive
  matrix[J, K] profiles;                         // 0/1 profiles, binary order
  int<lower=1> Q;                                // quadrature nodes
  vector[Q] nodes;                               // Gauss-Hermite nodes for N(0, 1)
  vector[Q] log_weights;                         // and log weights (summing to 1)
  array[J] int<lower=0> n1;                      // verified diseased
  array[J] int<lower=0> n0;                      // verified healthy
  array[J] int<lower=0> u;                       // unverified

  int prevalence_family;                         // prior on the prevalence
  vector[3] prevalence_params;
  // priors on the marginal probits m_dk: row 1 on the healthy, the prior of
  // 1 - Sp_k = Phi(m_1k); row 2 on the diseased, the prior of Se_k = Phi(m_2k)
  array[2, K] int accuracy_family;
  array[2, K] vector[3] accuracy_params;
  int loadings_family;                           // prior on the loadings b_dk
  vector[3] loadings_params;

  int<lower=0, upper=1> mnar;                    // 0 = MAR, 1 = MNAR
  int<lower=0> n_beta;                           // 0 (MAR or fixed beta), 1, or one per verified profile
  real beta_value;                               // fixed beta, used if mnar and n_beta == 0
  int alpha_family;                       // prior on alpha (used only if mnar)
  vector[3] alpha_params;
  int beta_family;                               // prior on beta (used only if n_beta > 0)
  vector[3] beta_params;
}

transformed data {
  array[J] int verified_pos = verified_positions(n1, n0);
  int J1 = max(verified_pos);
  if (n_beta > 1 && n_beta != J1) reject("n_beta must be 0, 1 or the number of verified profiles");
}

parameters {
  real logit_prevalence;
  // marginal probits: P(T_k = 1 | D = d) = Phi(m_dk). Sampling m instead of the
  // intercepts a = m sqrt(1 + b^2) removes the curved dependence between a and b
  array[2] vector[K] m;
  array[2] vector<lower=0>[positive_loadings ? K : 1] b_head;  // positive loadings
  array[2] vector[positive_loadings ? 0 : K - 1] b_tail;       // free loadings
  vector[mnar ? J1 : 0] alpha;
  vector[n_beta] beta;
}

transformed parameters {
  array[2] vector[K] b;                          // loadings
  array[2] vector[K] a;                          // intercepts
  for (d in 1:2) {
    b[d] = append_row(b_head[d], b_tail[d]);
    a[d] = m[d] .* sqrt(1 + square(b[d]));
  }
}

model {
  logit_prevalence ~ prior(prevalence_family, prevalence_params);
  for (d in 1:2) {
    for (k in 1:K) m[d][k] ~ prior(accuracy_family[d, k], accuracy_params[d, k]);
    // b_head is declared positive: Stan adds the Jacobian, the prior is truncated
    for (k in 1:K) b[d][k] ~ prior(loadings_family, loadings_params);
  }
  for (i in 1:num_elements(alpha)) alpha[i] ~ prior(alpha_family, alpha_params);
  for (i in 1:n_beta) beta[i] ~ prior(beta_family, beta_params);

  if (!prior_only) {
    target += cells_lp(logit_prevalence,
                       log_profiles(profiles, a[1], b[1], nodes, log_weights),
                       log_profiles(profiles, a[2], b[2], nodes, log_weights),
                       n1, n0, u, verified_pos, mnar, n_beta, beta_value, alpha, beta);
  }
}

generated quantities {
  real prevalence = inv_logit(logit_prevalence);
  vector[J] theta_healthy = exp(log_profiles(profiles, a[1], b[1], nodes, log_weights));
  vector[J] theta_diseased = exp(log_profiles(profiles, a[2], b[2], nodes, log_weights));
}
