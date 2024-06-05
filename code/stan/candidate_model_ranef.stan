// functions from brms 2.21.0
functions {
  /* hurdle lognormal log-PDF of a single response
   * Args:
   *   y: the response value
   *   mu: mean parameter of the lognormal distribution
   *   sigma: sd parameter of the lognormal distribution
   *   hu: hurdle probability
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real hurdle_lognormal_lpdf(real y, real mu, real sigma, real hu) {
    if (y == 0) {
      return bernoulli_lpmf(1 | hu);
    } else {
      return bernoulli_lpmf(0 | hu) +
             lognormal_lpdf(y | mu, sigma);
    }
  }
  /* hurdle lognormal log-PDF of a single response
   * logit parameterization of the hurdle part
   * Args:
   *   y: the response value
   *   mu: mean parameter of the lognormal distribution
   *   sigma: sd parameter of the lognormal distribution
   *   hu: linear predictor for the hurdle part
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real hurdle_lognormal_logit_lpdf(real y, real mu, real sigma, real hu) {
    if (y == 0) {
      return bernoulli_logit_lpmf(1 | hu);
    } else {
      return bernoulli_logit_lpmf(0 | hu) +
             lognormal_lpdf(y | mu, sigma);
    }
  }
  // hurdle lognormal log-CCDF and log-CDF functions
  real hurdle_lognormal_lccdf(real y, real mu, real sigma, real hu) {
    return bernoulli_lpmf(0 | hu) + lognormal_lccdf(y | mu, sigma);
  }
  real hurdle_lognormal_lcdf(real y, real mu, real sigma, real hu) {
    return log1m_exp(hurdle_lognormal_lccdf(y | mu, sigma, hu));
  }
  /* compute correlated group-level effects
  * Args:
  *   z: matrix of unscaled group-level effects
  *   SD: vector of standard deviation parameters
  *   L: cholesky factor correlation matrix
  * Returns:
  *   matrix of scaled group-level effects
  */
  matrix scale_r_cor(matrix z, vector SD, matrix L) {
    // r is stored in another dimension order than z
    return transpose(diag_pre_multiply(SD, L) * z);
  }
  /* softmax for a row vector
   * Args:
   *   y: (unconstrained) row vector
   * Returns:
   *   row vector simplex
   */
   row_vector softmax_row_vector(row_vector y) {
     return exp(y) / sum(exp(y));
   }
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  // split Y==0 and Y>0 for efficiency 
  int<lower=0> N_Ye0; // Number of y = 0 's
  array[N_Ye0] int indexes_Ye0; // Indexes of Y where y==0
  int<lower=0> N_Yg0; // Number of y > 0 's
  array[N_Yg0] int indexes_Yg0; // Indexes of Y where y > 0
  matrix[N, 1] X;  // design matrix
  matrix[N, 1] Xc;  // design matrix (centered)
  // data for group-level effects 
  int<lower=1> N_groups;  // number of grouping levels (sites)
  array[N] int<lower=1> J_group;  // grouping indicator per observation
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  // 1's and 0's to feed to bernoulli_lpmf; note y=0 --> bern(1|hu)
  array[N] int bern_01;
  bern_01[indexes_Ye0] = ones_int_array(N_Ye0);
  bern_01[indexes_Yg0] = zeros_int_array(N_Yg0);
  matrix[N_Yg0, 1] X_Yg0 = X[indexes_Yg0,];
  matrix[N_Yg0, 2] XInt_Yg0 = append_col(rep_vector(1.0, N_Yg0), X_Yg0);
  vector[N_Yg0] Y_Yg0 = Y[indexes_Yg0];
  array[N_Yg0] int<lower=1> J_group_Yg0 = J_group[indexes_Yg0];  // grouping indicator per observation
  vector[1] means_X;  // column means of X before centering
  means_X[1] = mean(X[,1]);
  matrix[N, 2] XcInt = append_col(rep_vector(1.0, N), Xc);  // centered version of X with an intercept
}
parameters {
  real b_b0;  // regression coefficient
  vector[1] b_IP;  // regression coefficient
  real<lower=0> sigma;  // dispersion parameter
  vector[1] b_hu;  // unscaled hurdle coefficients
  real Intercept_hu;  // temporary intercept for centered predictors
  vector<lower=0>[2] sd_grp_hu;  // group-level standard deviations
  matrix[2, N_groups] z_grp_hu;  // standardized group-level effects
  cholesky_factor_corr[2] L_grp_hu;  // cholesky factor of correlation matrix
}
transformed parameters {
  matrix[N_groups, 2] r_grp_hu; // actual group-level effects (unconstrained, diff from population effect)
  real lprior = 0;  // prior contributions to the log posterior
  // compute actual group-level effects
  r_grp_hu = scale_r_cor(z_grp_hu, sd_grp_hu, L_grp_hu);
  // priors
  lprior += normal_lpdf(b_b0 | 0, 1);
  lprior += normal_lpdf(b_IP | 0, 1);
  lprior += normal_lpdf(b_hu | 0, 1);
  lprior += student_t_lpdf(sigma | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += logistic_lpdf(Intercept_hu | 0, 1);
  lprior += student_t_lpdf(sd_grp_hu | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += lkj_corr_cholesky_lpdf(L_grp_hu | 1);
}
model {
  // likelihood including constants
  if (!prior_only) {
    vector[N_Yg0] mu;
    vector[N] hu;
    hu = Intercept_hu + Xc * b_hu;
    hu += rows_dot_product(r_grp_hu[J_group,], XcInt);
    mu = b_b0 + X_Yg0 * b_IP;
    target += bernoulli_logit_lpmf(bern_01 | hu);
    target += lognormal_lpdf(Y_Yg0 | mu, sigma);
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(to_vector(z_grp_hu));
}
generated quantities {
  // actual population-level intercept
  real b_hu_Intercept = Intercept_hu - dot_product(means_X, b_hu);
}

