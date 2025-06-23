library(ggplot2)
library(gridExtra)
library(MCuptake)
library(rstan)

country_code = "MWI"

## Configure Stan to use multiple cores. Fitting will take ~4 times as long if
## this is not done
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

fit_data = list(iso_code = country_code)
fit_data$pop_data = subset(mc_pop_data, ISO_Alpha_3==fit_data$iso_code & Year >= 1970 & Year <= 2025)
fit_data$svy_data = subset(mc_svy_data, ISO_Alpha_3==fit_data$iso_code)
fit_data$imis_fit = fit_mc_model(fit_data$pop_data, fit_data$svy_data)

# Save the model fit (e.g., "uga-fit.rds" for Uganda); this can be reloaded using fit_data = readRDS("uga-fit.rds")
save_rds = sprintf("%s-fit.rds", tolower(fit_data$iso_code))
saveRDS(fit_data, save_rds)

prev_tiff = sprintf("%s-mc-prev.tiff", tolower(fit_data$iso_code))
rate_tiff = sprintf("%s-mc-rate.tiff", tolower(fit_data$iso_code))
count_tiff = sprintf("%s-mc-count.tiff", tolower(fit_data$iso_code))
prev_csv = sprintf("%s-mc-prev.csv", tolower(fit_data$iso_code))

fit_data$pop_data = subset(mc_pop_data, ISO_Alpha_3==fit_data$iso_code & Year >= 1970 & Year <= 2050)
plot_fitted_mc_prev(prev_tiff, fit_data$imis_fit, fit_data$pop_data, fit_data$svy_data, 2000, 2050)
plot_fitted_mc_rates(rate_tiff, fit_data$imis_fit, fit_data$pop_data)
plot_fitted_mc_count(count_tiff, fit_data$imis_fit, fit_data$pop_data)
write_mc_prev(prev_csv, fit_data$imis_fit, fit_data$pop_data)
