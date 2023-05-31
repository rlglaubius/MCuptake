## Convert a parameter matrix (rows=samples, columns=parameters) to a list that
## can be used with (e.g.) mc_model_rate.
unpack_pars = function(par_matrix) {
  return(apply(par_matrix, 1, function(row_dat) {
    list(mc_uptake_1 = row_dat[1:5],
         mc_agedst_1 = row_dat[6:7],
         mc_uptake_2 = row_dat[8],
         mc_agedst_2 = row_dat[9:10])
  }))
}

## Male circumcision uptake rates over time
mc_model_rate = function(year, par) {
  z1 = 1.0 / (1.0 + exp(-par$mc_uptake_1[1] * (year - 1970 - par$mc_uptake_1[3])))
  z2 = 1.0 / (1.0 + exp(+par$mc_uptake_1[2] * (year - 1970 - par$mc_uptake_1[3])))
  mc_rate_1 = z1 * (2 * par$mc_uptake_1[4] * z2 + par$mc_uptake_1[5])
  mc_rate_2 = par$mc_uptake_2
  return(cbind(mc_rate_1, mc_rate_2))
}

## Age distribution of male circumcision
mc_model_dist = function(age, par) {
  return(cbind(dnbinom(age, size=par$mc_agedst_1[1], mu=par$mc_agedst_1[2]),
               dnbinom(age, size=par$mc_agedst_2[1], mu=par$mc_agedst_2[2])))
}

mc_model = function(year, age, par) {
  mc_rate = mc_model_rate(year, par)
  mc_dist = mc_model_dist(age,  par)
  return(outer(mc_rate[,1], mc_dist[,1]) + outer(mc_rate[,2], mc_dist[,2]))
}

model_sim = function(par, pop_data) {
  year = unique(pop_data$Year)
  ages = unique(pop_data$Age)
  num_a = length(ages)
  num_t = length(year)

  age_wide = reshape2::dcast(pop_data, Year~Age, value.var="Value")
  pop_sum = data.matrix(age_wide[,2:ncol(age_wide)])
  pop_unc = matrix(0, ncol=num_a, nrow=num_t) # uncircumcised population
  pop_crc = matrix(0, ncol=num_a, nrow=num_t) # circumcised population
  num_crc = matrix(0, ncol=num_a, nrow=num_t) # circumcisions performed

  mc_rate = mc_model(year, ages, par)
  mc_prob = 1.0 - exp(-mc_rate)

  ## Initialize the base-year population assuming that male circumcision
  ## prevalence has equilibrated according to base-year uptake rates. This
  ## assumes a stable population. If we have jitters in early MC prevalence, we
  ## could start simulation in 1950 instead of 1970 so that any early jitters
  ## are less likely to affect Spectrum inputs
  mc_init = rep(mc_prob[1,1], num_a)
  for (a in 2:num_a) {mc_init[a] = mc_init[a-1] + (1.0 - mc_init[a-1]) * mc_prob[1,a]}
  pop_crc[1,] = mc_init * pop_sum[1,]
  pop_unc[1,] = pop_sum[1,] - pop_crc[1,]
  num_crc[1,] = mc_init * pop_sum[1,]

  ## Initialize prevalence at age 0 ()
  pop_crc[,1] = pop_sum[,1] * mc_prob[,1]
  pop_unc[,1] = pop_sum[,1] - pop_crc[,1]

  ## Calculate MC prevalence from uptake and age 0 and base-year levels
  for (yi in 2:num_t) {
    sx = pop_sum[yi,2:num_a] / (pop_sum[yi-1,1:(num_a-1)] + 1e-16) # approximate survival from year-to-year populations. Padding to avoid divide-by-zero
    pop_crc[yi,2:num_a] = (pop_crc[yi-1,1:(num_a-1)] + mc_prob[yi,2:num_a] * pop_unc[yi-1,1:(num_a-1)]) * sx
    pop_unc[yi,2:num_a] = pop_sum[yi,2:num_a] - pop_crc[yi,2:num_a]
    num_crc[yi,2:num_a] = mc_prob[yi,2:num_a] * pop_unc[yi-1,1:(num_a-1)]
  }

  return(list(year = year, pop_sum = pop_sum, pop_unc = pop_unc, pop_crc = pop_crc, num_crc = num_crc))
}

sample_prior = function(n) {
  X = cbind(rexp(n, 0.5),           ## 1st mc rate, 1st slope parameter
            rexp(n, 0.5),           ## 1st mc rate, 2nd slope parameter
            rnorm(n, 2012-1970, 8), ## 1st mc rate, location parameter
            rexp(n, 0.5),           ## 1st mc rate, 1st height parameter
            rexp(n, 0.5),           ## 1st mc rate, 2nd height parameter

            rexp(n, 0.1),           ## 1st mc age distribution, negative binomial size
            runif(n, 0, 50),        ## 1st mc age distribution, negative binomial mean

            rexp(n, 8.0),

            rexp(n, 0.1),
            runif(n, 0, 50))
}

prior = function(X) {
  return(rowSums(
    cbind(
      dexp( X[,1], 0.5,  log=TRUE),
      dexp( X[,2], 0.5,  log=TRUE),
      dnorm(X[,3], 2012-1970, 8, log=TRUE),
      dexp( X[,4], 0.5,   log=TRUE),
      dexp( X[,5], 0.5,   log=TRUE),

      dexp( X[,6], 0.1,  log=TRUE),
      dunif(X[,7], 0, 50,log=TRUE),

      dexp( X[, 8], 8.0,   log=TRUE),

      dexp( X[,12], 0.1,  log=TRUE),
      dunif(X[,13], 0, 50,log=TRUE))))
}

likelihood = function(X, pop_data, svy_data) {
  P = prior(X)
  L = rep(-Inf, nrow(X))
  par_list = unpack_pars(X)
  for (k in which(is.finite(P))) {
    inputs = par_list[[k]]
    output = model_sim(inputs, pop_data)

    L[k] = sum(apply(svy_data[,c("age_min", "age_max", "year", "num_uncircumcised", "num_circumcised")], 1, function(obs) {
      a_ind = 1 + obs[1]:obs[2]
      y_ind = 1 + obs[3] - 1970
      numer = sum(output$pop_crc[y_ind,a_ind])
      denom = sum(output$pop_sum[y_ind,a_ind])
      prev  = numer / denom
      return(obs[5] * log(prev) + obs[4] * log(1.0 - prev))
    }))
  }
  return(L)
}

posterior = function(X) {
  return(prior(X) + likelihood(X))
}
