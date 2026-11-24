  // Prior densities for parameters declared on the unconstrained scale.
  //
  // Every scalar parameter of the baysesp models is declared on the real line:
  // a probability p is represented by x = logit(p). The prior is coded in R by
  // prior_code() as a family code and a vector of three parameters:
  //   1 = normal(mean, sd)                 params = [mean, sd, 0]
  //   2 = student_t(df, location, scale)   params = [df, location, scale]
  //   3 = beta(shape1, shape2)             params = [shape1, shape2, 0]
  // Families 1 and 2 are densities of x itself (for a probability: of its
  // logit). Family 3 is a density of p = inv_logit(x), so the log Jacobian of
  // the transformation, log p + log(1 - p), is added.

  real prior_lpdf(real x, int family, vector params) {
    if (family == 1) {
      return normal_lpdf(x | params[1], params[2]);
    } else if (family == 2) {
      return student_t_lpdf(x | params[1], params[2], params[3]);
    } else if (family == 3) {
      return beta_lpdf(inv_logit(x) | params[1], params[2])
             + log_inv_logit(x) + log1m_inv_logit(x);
    }
    reject("prior_lpdf: unknown prior family ", family);
  }

  // Same prior for every element of a vector.
  real prior_vector_lpdf(vector x, int family, vector params) {
    real lp = 0;
    for (i in 1:num_elements(x)) {
      lp += prior_lpdf(x[i] | family, params);
    }
    return lp;
  }
