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
  real<lower=0> r1_slope1;
  real<lower=0> r1_slope2;
  real<lower=0> r1_center;
  real<lower=0> r1_theta1;
  real<lower=0> r1_theta2;

  // Age distribution used with the first rate
  real<lower=0> a1_size;
  real<lower=0,upper=50> a1_mean;

  real<lower=0> r2_rate;

  // Age distribution used with the second rate
  real<lower=0> a2_size;
  real<lower=0,upper=50> a2_mean;
}

transformed parameters {
  real z1;
  real z2;
  real<lower=0> rate1[n_pop_yrs];
  real<lower=0> rate2[n_pop_yrs];
  real<lower=0> dist1[n_pop_age];
  real<lower=0> dist2[n_pop_age];

  real<lower=0,upper=1> uptake_prop[n_pop_yrs, n_pop_age];
  real<lower=0,upper=1> mc_prev[n_pop_yrs,n_pop_age];

  real<lower=0, upper=1> svy_mc_prev[n_svy_obs];
  real<lower=0> numer;
  real<lower=0> denom;

  // Calculate time trends in rates
  for (t in 1:n_pop_yrs) {
    z1 = 1.0 / (1.0 + exp(-r1_slope1 * (pop_year[t] - 1970 - r1_center)));
    z2 = 1.0 / (1.0 + exp(+r1_slope2 * (pop_year[t] - 1970 - r1_center)));
    rate1[t] = z1 * (2 * r1_theta1 * z2 + r1_theta2);
    rate2[t] = r2_rate;
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
      numer += pop_sum[svy_year[k]-1970+1,a+1] * mc_prev[svy_year[k]-1970+1,a+1];
      denom += pop_sum[svy_year[k]-1970+1,a+1];
    }
    svy_mc_prev[k] = numer / denom;
  }
}

model {
  // Prior
  r1_slope1 ~ exponential(0.5);
  r1_slope2 ~ exponential(0.5);
  r1_center ~ normal(2012-1970, 4.0);
  r1_theta1 ~ exponential(2.0);
  r1_theta2 ~ exponential(0.5);

  a1_size ~ exponential(0.1);
  a1_mean ~ uniform(0, 50);

  r2_rate ~ exponential(8.0);

  a2_size ~ exponential(0.1);
  a2_mean ~ uniform(0, 50);

  // Likelihood
  for (k in 1:n_svy_obs) {
    target += svy_n_mc[k] * log(svy_mc_prev[k]) + svy_n_no[k] * log(1.0 - svy_mc_prev[k]);
  }
}

