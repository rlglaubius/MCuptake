#' Run a model fit for one country
#'
#' Estimate male circumcision prevalence and produce a diagnostic plot showing
#' goodness-of-fit to available survey data, as well as point estimates of
#' male circumcision prevalence trends by five-year age group (15-19, 20-24, ..., 45-49)
#' @param isocode ISO-3166 3-letter alphabetic code for the country to be fitted (e.g., "UGA" for Uganda)
#' @param tiffname Figure name for output plot showing goodness of fit. This should end in .tiff
#' @param csvname File name for an output CSV of male circumcision prevalence point estimates
#' @export
estimate_mc_prev = function(isocode, tiffname, csvname) {
  pop_data = subset(mc_pop_data, ISO_Alpha_3==iso_code)
  if (nrow(pop_data) == 0) {
    error(sprintf("No population data found for %s", isocode))
  }

  svy_data = subset(mc_svy_data, ISO_Alpha_3==iso_code)
  if (nrow(svy_data) == 0) {
    error(sprintf("No survey data found for %s", isocode))
  }

  imis_fit = fit_model(pop_data, svy_data)
  plot_fitted_prev(tiffname, imis_fit, pop_data, svy_data)
  write_mc_prev(csvname, imis_fit, pop_data)
}

#' Fit the male circumcision uptake model for one country
#' @param pop_data A long data frame of male population sizes by year and single age for the selected country
#' @param svy_data A data frame of male circumcision prevalence estimates from household surveys for the selected country
#' @return An IMIS model fit object (See details section)
#' @section Details:
#'
#' The IMIS model fit object is a list of the following:
#' \enumerate{
#' \item{resample} A matrix of posterior parameter values, one row per sample, one column per parameter
#' \item{center} IMIS constructs a gaussian mixture model approximation to the posterior. "center" is a matrix of
#' modes of the gaussian components of that mixture model
#' \item{stat} Convergence monitoring statistics from IMIS iterations. These are also displayed in the console
#' when fit_model is running
#' \item{prior} log prior density at each resample
#' \item{lhood} log likelihood at each resample
#' }
#'
#' @export
fit_model = function(pop_data, svy_data) {
  lhood = function(X) {likelihood(X, pop_data, svy_data)}
  imis_fit = log_imis(prior, lhood, sample_prior, 500, 1000, 1000) # list(resample=sample_prior(100)
  imis_fit$prior = prior(imis_fit$resample)
  imis_fit$lhood = lhood(imis_fit$resample)
  return(imis_fit)
}
