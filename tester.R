library(MCuptake)

iso_code = "UGA"
pop_data = subset(mc_pop_data, ISO_Alpha_3==iso_code & Year >= 1970 & Year <= 2025)
svy_data = subset(mc_svy_data, ISO_Alpha_3==iso_code)

imis_fit = fit_mc_model(pop_data, svy_data)
plot_fitted_mc_prev("fig-mc-prev-uga.tiff", imis_fit, pop_data, svy_data)
write_mc_prev("csv-mc-prev-uga.csv", imis_fit, pop_data)

# estimate_mc_prev("UGA", "fig-mc-prev-uga.tiff", "csv-mc-prev-uga.csv")
