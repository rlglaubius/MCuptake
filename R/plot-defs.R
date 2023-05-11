#' Plot fitted male circumcision prevalence
#'
#' Plot fitted male circumcision prevalence for a single country as a TIFF file
#' @param tiffname Output figure filename (tiff format)
#' @param imis_fit IMIS model fit object
#' @param pop_data A long data frame of male population sizes by year and single
#'   age for the selected country
#' @param svy_data A data frame of male circumcision prevalence estimates from
#'   household surveys for the selected country
#' @export
plot_fitted_prev = function(tiffname, imis_fit, pop_data, svy_data) {
  par_list = apply(imis_fit$resample, 1, function(row_dat) {
    list(mc_uptake_1 = row_dat[ 1:6 ],
         mc_agedst_1 = row_dat[ 7:8 ],
         mc_uptake_2 = row_dat[ 9:12],
         mc_agedst_2 = row_dat[13:14])
  })

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
    xlim(c(2000,2025)) + ylim(c(0,100)) +
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
          axis.text.x = element_text(color="#000000", angle=45, hjust=1))
  ggsave(tiffname, compression="lzw", dpi=600, width=2*3.42, height=2*2.44)
}

#' Write male circumcision prevalence trends to a CSV file
#' @param csvname File name for an output CSV of male circumcision prevalence point estimates
#' @param imis_fit IMIS model fit object
#' @param pop_data A long data frame of male population sizes by year and single
#'   age for the selected country
#' @export
write_mc_prev = function(csvname, imis_fit, pop_data) {
  ind_best = which.max(imis_fit$prior + imis_fit$lhood)
  par_best = list(mc_uptake_1 = imis_fit$resample[ind_best, 1:6 ],
                  mc_agedst_1 = imis_fit$resample[ind_best, 7:8 ],
                  mc_uptake_2 = imis_fit$resample[ind_best, 9:12],
                  mc_agedst_2 = imis_fit$resample[ind_best,13:14])
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

