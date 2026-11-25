  // Log-likelihood of the aggregated table for the structured models
  // (specification, section 3.4), from the prevalence on the logit scale and
  // the log distributions of the profiles in the healthy and in the diseased.
  // Verification: ignorable under MAR (mnar = 0); under MNAR
  // logit(lambda_{t,d}) = alpha_t + beta d, with beta fixed (n_beta = 0), one
  // beta (n_beta = 1) or one per verified profile. Never-verified profiles
  // (verified_pos = 0) have lambda = 0.

  real cells_lp(real logit_prevalence, vector log_theta0, vector log_theta1,
                array[] int n1, array[] int n0, array[] int u,
                array[] int verified_pos, int mnar, int n_beta, real beta_value,
                vector alpha, vector beta) {
    real log_p1 = log_inv_logit(logit_prevalence);
    real log_p0 = log1m_inv_logit(logit_prevalence);
    real lp = 0;
    for (j in 1:num_elements(n1)) {
      real lp1 = log_p1 + log_theta1[j];         // log P(D = 1, T = t)
      real lp0 = log_p0 + log_theta0[j];         // log P(D = 0, T = t)
      int a = verified_pos[j];
      if (a == 0) {
        lp += u[j] * log_sum_exp(lp1, lp0);
      } else if (mnar == 0) {
        lp += n1[j] * lp1 + n0[j] * lp0 + u[j] * log_sum_exp(lp1, lp0);
      } else {
        real b = n_beta == 0 ? beta_value : beta[n_beta == 1 ? 1 : a];
        real v1 = alpha[a] + b;                  // logit lambda_{t,1}
        real v0 = alpha[a];                      // logit lambda_{t,0}
        lp += n1[j] * (lp1 + log_inv_logit(v1)) + n0[j] * (lp0 + log_inv_logit(v0));
        lp += u[j] * log_sum_exp(lp1 + log1m_inv_logit(v1), lp0 + log1m_inv_logit(v0));
      }
    }
    return lp;
  }

  // Position of each profile among those with verified subjects, 0 if never verified.
  array[] int verified_positions(array[] int n1, array[] int n0) {
    int J = num_elements(n1);
    array[J] int pos;
    int a = 0;
    for (j in 1:J) {
      if (n1[j] + n0[j] > 0) {
        a += 1;
        pos[j] = a;
      } else {
        pos[j] = 0;
      }
    }
    return pos;
  }
