// baysesp: modello saturo in parametrizzazione trasparente.
// Parametri identificati: omega (frequenze dei profili), v (quote di verificati),
// phi1 (proporzione di malati tra i verificati di ciascun profilo).
// Parametri di sensibilita', non informati dai dati e con priori esplicite:
// q (proporzione di malati nei profili mai verificati) e beta (MNAR).

data {
  int<lower=1> K;
  int<lower=2> J;
  matrix[K, J] A;                          // matrice dei profili
  int<lower=0> R;
  matrix[R, J] Rmat;                       // vettori r_R disposti per riga

  array[J] int<lower=0> n1;                // verificati malati
  array[J] int<lower=0> n0;                // verificati sani
  array[J] int<lower=0> u;                 // non verificati
  array[J] int<lower=0, upper=1> mai_ver;  // 1 se il profilo non e' mai verificato

  int<lower=0, upper=1> mnar;              // 0 = MAR, 1 = MNAR con parametro beta
  vector<lower=0>[J] a_q;                  // priori Beta(a_q, b_q) su q_t
  vector<lower=0>[J] b_q;                  //   (usate solo per i profili mai verificati)
  real beta_media;
  real<lower=0> beta_sd;
}

transformed data {
  int J0 = sum(mai_ver);                   // numero di profili mai verificati
  int J1 = J - J0;
  array[J1] int idx_v;                     // indici dei profili verificati
  array[J0] int idx_0;                     // indici dei profili mai verificati
  array[J] int Nt;
  {
    int a = 1;
    int b = 1;
    for (j in 1:J) {
      Nt[j] = n1[j] + n0[j] + u[j];
      if (mai_ver[j] == 1) { idx_0[b] = j; b += 1; }
      else { idx_v[a] = j; a += 1; }
    }
  }
}

parameters {
  simplex[J] omega;                        // P(T = t)
  vector<lower=0, upper=1>[J1] v;          // P(V = 1 | T = t), profili verificati
  vector<lower=0, upper=1>[J1] phi1;       // P(D = 1 | T = t, V = 1)
  vector<lower=0, upper=1>[J0] q;          // P(D = 1 | T = t), profili mai verificati
  vector[mnar] beta;                       // log-odds ratio di verifica malati/sani
}

model {
  // priori sui parametri identificati: deboli, la verosimiglianza domina
  omega ~ dirichlet(rep_vector(1, J));
  v ~ beta(1, 1);
  phi1 ~ beta(1, 1);
  // priori sui parametri di sensibilita': determinano da sole le loro posteriori
  for (b in 1:J0) q[b] ~ beta(a_q[idx_0[b]], b_q[idx_0[b]]);
  if (mnar) beta ~ normal(beta_media, beta_sd);

  // verosimiglianza: multinomiale sui profili e binomiali indipendenti
  Nt ~ multinomial(omega);
  for (a in 1:J1) {
    int j = idx_v[a];
    n1[j] + n0[j] ~ binomial(Nt[j], v[a]);
    n1[j] ~ binomial(n1[j] + n0[j], phi1[a]);
  }
}

generated quantities {
  vector[J] psi;                           // P(D = 1 | T = t)
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
