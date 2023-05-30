data {
  int n_pop_yrs;
  int n_pop_age;
  int<lower=0>  pop_year[n_pop_yrs];
  real<lower=0> pop_sum[n_pop_yrs, n_pop_age];

  int n_svy_obs;
  int<lower=0> svy_age_min[n_svy_obs];
  int<lower=0> svy_age_max[n_svy_obs];
  int<lower=0> svy_year[n_svy_obs];
  real<lower=0> svy_n_mc[n_svy_obs]; // number circumcised
  real<lower=0> svy_n_no[n_svy_obs]; // number uncircumcised
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  // The first rate (r1) is defined as a probit-transformed logistic trend. The probit transform allows the
  // rate to peak and decline
  real<upper=0> r1_llim; // left limit
  real<lower=0> r1_rlim; // right limit
  real<lower=0> r1_slope;
  real r1_midpt;
  real<lower=0> r1_max; // rate maximum
  real<lower=0> r1_std; // probit standard deviation

  // Age distribution used with the first rate
  real<lower=0> a1_size;
  real<lower=0,upper=50> a1_mean;

  // The second rate (r2) is defined on its natural scale, so left and right limits must be non-negative
  real<lower=0> r2_llim; // left limit
  real<lower=0> r2_rlim; // right limit
  real<lower=0> r2_slope;
  real r2_midpt;

  // Age distribution used with the second rate
  real<lower=0> a2_size;
  real<lower=0,upper=50> a2_mean;
}

transformed parameters {
  real z1;
  real<lower=0> z2;
  real<lower=0> rate1[n_pop_yrs];
  real<lower=0> rate2[n_pop_yrs];
  real<lower=0> dist1[n_pop_age];
  real<lower=0> dist2[n_pop_age];

  real<lower=0,upper=1> uptake_prop[n_pop_yrs, n_pop_age];
  real<lower=0,upper=1> mc_prev[n_pop_yrs,n_pop_age];

  real<lower=0, upper=1> svy_mc_prev[n_svy_obs];
  real<lower=0> numer;
  real<lower=0> denom;

  // real<lower=0,upper=1> mc_init[n_pop_age];
  // real<lower=0> pop_crc[n_pop_yrs, n_pop_age];
  // real<lower=0> pop_unc[n_pop_yrs, n_pop_age];

  // Calculate time trends in rates
  for (t in 1:n_pop_yrs) {
    z1 = r1_llim + (r1_rlim - r1_llim) / (1.0 + exp(-r1_slope * (pop_year[t] - 1970 - r1_midpt)));
    z2 = r2_llim + (r2_rlim - r2_llim) / (1.0 + exp(-r2_slope * (pop_year[t] - 1970 - r2_midpt)));
    rate1[t] = r1_max * exp(-0.5 * (z1 / r1_std) * (z1 / r1_std));
    rate2[t] = z2;
  }

  // Calculate rate age distributions
  for (a in 1:n_pop_age) {
    dist1[a] = exp(neg_binomial_2_lpmf(a-1 | a1_mean, a1_size));
    dist2[a] = exp(neg_binomial_2_lpmf(a-1 | a2_mean, a2_size));
  }

  for (t in 1:n_pop_yrs) {
    for (a in 1:n_pop_age) {
      uptake_prop[t,a] = 1.0 - exp(-rate1[t] * dist1[a] - rate2[t] * dist2[a]);
    }
  }

  // Initialize age 0
  for (t in 1:n_pop_yrs) {
    mc_prev[t,1] = uptake_prop[t,1];
  }

  // Initialize the first year
  for (a in 2:n_pop_age) {
    mc_prev[1,a] = mc_prev[1,a-1] + (1.0 - mc_prev[1,a-1]) * uptake_prop[1,a];
  }

  // Initialize interior
  for (t in 2:n_pop_yrs) {
    for (a in 2:n_pop_age) {
      mc_prev[t,a] = mc_prev[t-1,a-1] + (1.0 - mc_prev[t-1,a-1]) * uptake_prop[t,a];
    }
  }

  // Calculate model estimates for survey data points
  for (k in 1:n_svy_obs) {
    numer = 0.0;
    denom = 0.0;
    for (a in svy_age_min[k]:svy_age_max[k]) {
      numer += pop_sum[svy_year[k]-1970+1,a] * mc_prev[svy_year[k]-1970+1,a];
      denom += pop_sum[svy_year[k]-1970+1,a];
    }
    svy_mc_prev[k] = numer / denom;
  }

  // The version below tries to take population dynamics into account. This may not be
  // necessary, and seems to cause numerical problems.
  // // Initialize the first year
  // mc_init[1] = uptake_prop[1,1];
  // for (a in 2:n_pop_age) {
  //   mc_init[a] = mc_init[a-1] + (1.0 - mc_init[a-1]) * uptake_prop[1,a];
  //   pop_crc[1,a] = mc_init[a] * pop_sum[1,a];
  //   pop_unc[1,a] = pop_sum[1,a] - pop_crc[1,a];
  // }
  //
  // // Initialize age 0
  // for (t in 2:n_pop_yrs) {
  //   pop_crc[t,1] = pop_sum[t,1] * uptake_prop[t,1];
  //   pop_unc[t,1] = pop_sum[t,1] - pop_crc[t,1];
  // }
  //
  // for (t in 2:n_pop_yrs) {
  //   for (a in 2:n_pop_age) {
  //     pop_crc[t,a] = (pop_crc[t-1,a-1] + pop_unc[t-1,a-1] * uptake_prop[t,a]) * (pop_sum[t,a] / pop_sum[t-1,a-1]);
  //     pop_unc[t,a] = pop_sum[t,a] - pop_crc[t,a];
  //   }
  // }
}

model {
  // Prior
  r1_llim ~ normal(0,1);
  r1_rlim ~ normal(0,1);
  r1_slope ~ exponential(0.5);
  r1_midpt ~ normal(2012-1970, 8.0);
  r1_max ~ exponential(1.0);
  r1_std ~ exponential(1.0);

  a1_size ~ exponential(0.1);
  a1_mean ~ uniform(0, 50);

  r2_llim ~ exponential(8.0);
  r2_rlim ~ exponential(8.0);
  r2_slope ~ exponential(0.5);
  r2_midpt ~ normal(2012-1970, 8.0);

  a2_size ~ exponential(0.1);
  a2_mean ~ uniform(0, 50);

  // Likelihood
  for (k in 1:n_svy_obs) {
    target += svy_n_mc[k] * log(svy_mc_prev[k]) + svy_n_no[k] * log(1.0 - svy_mc_prev[k]);
  }
}

