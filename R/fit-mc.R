#' Run a model fit for one country
#'
#' Estimate male circumcision prevalence and produce a diagnostic plot showing
#' goodness-of-fit to available survey data, as well as point estimates of
#' male circumcision prevalence trends by five-year age group (15-19, 20-24, ..., 45-49)
#' @param isocode ISO-3166 3-letter alphabetic code for the country to be fitted (e.g., "UGA" for Uganda)
#' @param tiffname Figure name for output plot showing goodness of fit. This should end in .tiff
#' @param csvname File name for an output CSV of male circumcision prevalence point estimates
#' @param first_year First year of estimates
#' @param final_year Final year of estimates
#' @export
estimate_mc_prev = function(isocode, tiffname, csvname, year_first=1970, year_final=2030) {
  if (year_first < 1950 | year_first > 1970) {
    error("year_first cannot be before 1950 or after 1970")
  }

  if (year_final < 2025 | year_final > 2100) {
    error("year_final cannot be before 2025 or after 2100")
  }

  pop_data = subset(mc_pop_data, ISO_Alpha_3==iso_code & Year >= year_first & Year <= year_final)
  if (nrow(pop_data) == 0) {
    error(sprintf("No population data found for %s", isocode))
  }

  svy_data = subset(mc_svy_data, ISO_Alpha_3==iso_code)
  if (nrow(svy_data) == 0) {
    error(sprintf("No survey data found for %s", isocode))
  }

  imis_fit = fit_mc_model(pop_data, svy_data)
  plot_fitted_mc_prev(tiffname, imis_fit, pop_data, svy_data)
  write_mc_prev(csvname, imis_fit, pop_data)
}

#' Fit the male circumcision uptake model for one country
#' @param pop_data A long data frame of male population sizes by year and single age for the selected country
#' @param svy_data A data frame of male circumcision prevalence estimates from household surveys for the selected country
#' @return An model fit object (See details section)
#' @section Details:
#'
#' The model fit object is a list of the following:
#' \enumerate{
#' \item{resample} A matrix of posterior parameter values, one row per sample, one column per parameter
#' \item{posterior} Posterior density at each sample. This is calculated by Stan, and accounts for
#' non-negativity constraints on parameters, which will cause the values to differ from those calculated
#' using \code{prior()} + \code{likelihood()}.
#' \item{raw} Raw output from \code{stan()}. Useful for checking convergence diagnostics.
#' }
#'
#' @export
fit_mc_model = function(pop_data, svy_data) {
  stan_data = list(
    n_pop_yrs = length(unique(pop_data$Year)),
    n_pop_age = length(unique(pop_data$Age)),
    pop_year = unique(pop_data$Year),
    pop_sum = pop_wide,

    n_svy_obs = nrow(svy_data),
    svy_age_min = svy_data$age_min,
    svy_age_max = svy_data$age_max,
    svy_year = svy_data$year,
    svy_n_mc = svy_data$num_circumcised,
    svy_n_no = svy_data$num_uncircumcised)

  ## Stan produces many more outputs. We retain just the MC uptake rate
  ## parameters
  keep_par = c("r1_slope1", "r1_slope2", "r1_center", "r1_theta1", "r1_theta2",
               "a1_size", "a1_mean",
               "r2_llim", "r2_rlim", "r2_slope", "r2_midpt",
               "a2_size", "a2_mean")

  stan_fit = stan(
    file = "R/mc-model.stan",
    data = stan_data,
    iter = 2000,
    chains = 4,
    include = TRUE,
    pars = keep_par)

  post_samp = as.data.frame(stan_fit)

  ## We repackage the Stan outputs for compatibility of the IMIS-based code
  ## version
  fit_obj = list()
  fit_obj$resample = as.matrix(post_samp[,keep_par]) # drop lp__
  fit_obj$posterior = post_samp$lp__
  fit_obj$raw = stan_fit
  return(fit_obj)
}
