library(gridExtra)
library(MCuptake)

iso_code = "UGA"
pop_data = subset(mc_pop_data, ISO_Alpha_3==iso_code & Year >= 1970 & Year <= 2025)
svy_data = subset(mc_svy_data, ISO_Alpha_3==iso_code)

prev_tiff = sprintf("%s-mc-prev.tiff", tolower(iso_code))
rate_tiff = sprintf("%s-mc-rate.tiff", tolower(iso_code))
prev_csv = sprintf("%s-mc-prev.csv", tolower(iso_code))

imis_fit = fit_mc_model(pop_data, svy_data)
plot_fitted_mc_prev(prev_tiff, imis_fit, pop_data, svy_data)
plot_fitted_mc_rates(rate_tiff, imis_fit, pop_data)
write_mc_prev(prev_csv, imis_fit, pop_data)
