
#' Plot fitted male circumcision prevalence
#'
#' Plot fitted male circumcision prevalence for a single country as a TIFF file
#' @param tiffname Output figure filename (tiff format)
#' @param imis_fit IMIS model fit object
#' @param pop_data A long data frame of male population sizes by year and single
#'   age for the selected country
#' @param svy_data A data frame of male circumcision prevalence estimates from
#'   household surveys for the selected country
#' @param year_first First year to plot
#' @param year_final Final year to plot
#' @export
plot_fitted_mc_prev = function(tiffname, imis_fit, pop_data, svy_data, year_first=2000, year_final=2025) {
  par_list = unpack_pars(imis_fit$resample)
  mod_list = lapply(par_list, function(par) {model_sim(par, pop_data)})
  ind_best = which.max(imis_fit$prior + imis_fit$lhood)
  years = mod_list[[1]]$year

  ## Calculate MC prevalence for each parameter set and survey datapoint
  svy_data$AgeGroup = sprintf("%d-%d", svy_data$age_min, svy_data$age_max)
  svy_summ = svy_data[,c("year", "AgeGroup", "prop_circumcised", "ci_lower", "ci_upper")]
  colnames(svy_summ) = c("Year", "AgeGroup", "Value", "Lower", "Upper")

  mod_summ = plyr::ddply(svy_data, .variables=c("AgeGroup"), function(df) {
    age_min = df$age_min[1]
    age_max = df$age_max[1]

    mc_prev = sapply(mod_list, function(mod) {
      numer = rowSums(mod$pop_crc[,age_min:age_max + 1])
      denom = rowSums(mod$pop_sum[,age_min:age_max + 1])
      return(numer / denom)
    })

    bounds = apply(mc_prev, 1, function(row_dat) {quantile(row_dat, c(0.025, 0.975))})
    data.frame(AgeGroup = df$AgeGroup[1],
               Year     = years,
               Value    = mc_prev[,ind_best],
               Lower    = bounds[1,],
               Upper    = bounds[2,])
  })

  plot_data = dplyr::bind_rows(list(Model = mod_summ, Data = svy_summ), .id="Source")

  ggplot(plot_data, aes(x=Year, y=100*Value, ymin=100*Lower, ymax=100*Upper, color=Source, fill=Source)) +
    geom_ribbon(data=plot_data[plot_data$Source=="Model",], show.legend=FALSE, alpha=0.2, color=NA) +
    geom_line(data=plot_data[plot_data$Source=="Model",]) +
    geom_point(data=plot_data[plot_data$Source=="Data",]) +
    geom_pointrange(data=plot_data[plot_data$Source=="Data",], size=0, show.legend=FALSE) +
    xlim(c(year_first,year_final)) + ylim(c(0,100)) +
    facet_wrap(~AgeGroup, nrow=2) +
    ylab("Male circumcision prevalence, %") +
    theme_bw() +
    theme(legend.position="right",
          legend.margin = margin(t=0, b=0, l=0, r=0),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(color="#a0a0a0", linewidth=0.05),
          plot.margin = margin(t=0, b=0, l=0.05, r=0.05, unit="cm"),
          strip.background = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(color="#000000", hjust=1, angle=45),
          axis.text.y = element_text(color="#000000"))
  ggsave(tiffname, compression="lzw", dpi=600, width=2*3.42, height=2*2.44)
}

#' Plot fitted male circumcision uptake rates
#'
#' Plot fitted male circumcision uptake rates for a single country as a TIFF file
#' @param tiffname Output figure filename (tiff format)
#' @param imis_fit IMIS model fit object
#' @param pop_data A long data frame of male population sizes by year and single
#'   age for the selected country
#' @export
plot_fitted_mc_rates = function(tiffname, imis_fit, pop_data) {
  matrix_row_summary = function(m, best) {
    data.frame(Value = m[,best],
               Lower = apply(m, 1, function(rdat) {quantile(rdat, c(0.025))}),
               Upper = apply(m, 1, function(rdat) {quantile(rdat, c(0.975))}))
  }

  par_list = unpack_pars(imis_fit$resample)
  ind_best = which.max(imis_fit$prior + imis_fit$lhood)

  age = 0:100
  yrs = unique(pop_data$Year)
  rate_list = lapply(par_list, function(par) {mc_model_rate(yrs, par)})
  dist_list = lapply(par_list, function(par) {mc_model_dist(age, par)})

  rate_1_mtrx = sapply(rate_list, function(r) {r[,1]})
  rate_2_mtrx = sapply(rate_list, function(r) {r[,2]})
  dist_1_mtrx = sapply(dist_list, function(d) {d[,1]})
  dist_2_mtrx = sapply(dist_list, function(d) {d[,2]})

  rate_1_summ = matrix_row_summary(rate_1_mtrx, ind_best)
  rate_2_summ = matrix_row_summary(rate_2_mtrx, ind_best)
  dist_1_summ = matrix_row_summary(dist_1_mtrx, ind_best)
  dist_2_summ = matrix_row_summary(dist_2_mtrx, ind_best)

  rate_1_summ$Year = yrs
  rate_2_summ$Year = yrs
  dist_1_summ$Age = age
  dist_2_summ$Age = age

  plot_rate = dplyr::bind_rows(list(Rate1 = rate_1_summ, Rate2 = rate_2_summ), .id="Component")
  plot_dist = dplyr::bind_rows(list(Rate1 = dist_1_summ, Rate2 = dist_2_summ), .id="Component")

  pr = ggplot(plot_rate, aes(x=Year, y=Value, ymin=Lower, ymax=Upper, color=Component, fill=Component)) +
    geom_ribbon(alpha=0.2, color=NA) +
    geom_line() +
    theme_bw() +
    ylab("MC uptake trend") +
    theme(legend.position="top",
          legend.margin = margin(t=0, b=0, l=0, r=0),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(color="#a0a0a0", linewidth=0.05),
          plot.margin = margin(t=0, b=0, l=0.05, r=0.1, unit="cm"),
          strip.background = element_blank(),
          axis.text = element_text(color="#000000"))

  pd = ggplot(plot_dist, aes(x=Age, y=Value, ymin=Lower, ymax=Upper, color=Component, fill=Component)) +
    geom_ribbon(alpha=0.2, color=NA) +
    geom_line() +
    theme_bw() +
    ylab("MC uptake age distribution") +
    theme(legend.position="top",
          legend.margin = margin(t=0, b=0, l=0, r=0),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(color="#a0a0a0", linewidth=0.05),
          plot.margin = margin(t=0, b=0, l=0.1, r=0.05, unit="cm"),
          strip.background = element_blank(),
          axis.text = element_text(color="#000000"))

  ggsave(tiffname, plot=grid.arrange(pr, pd, ncol=2, nrow=1), compression="lzw", dpi=600, width=2*3.42, height=2.44)
}

#' Plot numbers of male circumcisions done
#'
#' Plot numbers of male circumcisions done for a single country as a TIFF file
#' @param tiffname Output figure filename (tiff format)
#' @param imis_fit IMIS model fit object
#' @param pop_data A long data frame of male population sizes by year and single
#'   age for the selected country
#' @param year_first First year to plot
#' @param year_final Final year to plot
#' @export
plot_fitted_mc_count = function(tiffname, imis_fit, pop_data, svy_data, year_first=2000, year_final=2025) {
  par_list = unpack_pars(imis_fit$resample)
  ind_best = which.max(imis_fit$prior + imis_fit$lhood)
  year = unique(pop_data$Year)
  ages = unique(pop_data$Age)

  mc_count = sapply(par_list, function(par) {
    vals = model_sim(par, pop_data)
    return(rowSums(vals$num_crc))
  })

  bounds = apply(mc_count, 1, function(row_dat) {quantile(row_dat, c(0.025, 0.975))})
  plot_data = data.frame(Year = year,
                         Value = mc_count[,ind_best],
                         Lower = bounds[1,],
                         Upper = bounds[2,])

  ggplot(plot_data, aes(x=Year, y=1e-3*Value, ymin=1e-3*Lower, ymax=1e-3*Upper)) +
    geom_ribbon(alpha=0.2, color=NA) +
    geom_line() +
    xlim(c(year_first,year_final)) +
    ylab("Male circumcisions, 1000s") +
    expand_limits(y=0) +
    theme_bw()
  ggsave(tiffname, compression="lzw", dpi=600, width=3.42, height=2.44)
}

#' Write male circumcision prevalence trends to a CSV file
#' @param csvname File name for an output CSV of male circumcision prevalence point estimates
#' @param imis_fit IMIS model fit object
#' @param pop_data A long data frame of male population sizes by year and single
#'   age for the selected country
#' @export
write_mc_prev = function(csvname, imis_fit, pop_data) {
  par_list = unpack_pars(imis_fit$resample)

  ind_best = which.max(imis_fit$prior + imis_fit$lhood)
  par_best = par_list[[ind_best]]
  mod_best = model_sim(par_best, pop_data)

  age_groups = data.frame(age_min  = seq(15, 45, 5),
                          age_max  = seq(19, 49, 5))
  age_groups$AgeGroup = sprintf("%d-%d", age_groups$age_min, age_groups$age_max)
  mod_prev = plyr::ddply(age_groups, .variables=c("AgeGroup"), function(df) {
    age_min = df$age_min[1]
    age_max = df$age_max[1]
    numer = rowSums(mod_best$pop_crc[,age_min:age_max + 1])
    denom = rowSums(mod_best$pop_sum[,age_min:age_max + 1])
    data.frame(Age   = df$AgeGroup[1],
               Year  = mod_best$year,
               Value = numer / denom)
  })
  mod_prev$AgeGroup = NULL
  mod_wide = reshape2::dcast(mod_prev, Age~Year, value.var = "Value")
  write.csv(mod_wide, csvname, row.names=FALSE)
}

