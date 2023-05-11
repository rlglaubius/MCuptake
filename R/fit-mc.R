source("R/model-defs.R") # TODO: Remove
source("R/logimis.R")    # TODO: Remove
library(ggplot2)

#' Fit the male circumcision uptake model
#' @param pop_data A long data frame of male population sizes by year and single age
#' @param svy_data A data frame of male circumcision prevalence estimates from household surveys
#' @export
fit_model = function(pop_data, svy_data) {
  lhood = function(X) {likelihood(X, pop_data, svy_data)}
  imis_fit = log_imis(prior, lhood, sample_prior, 500, 1000, 1000)
  imis_fit$prior = prior(imis_fit$resample)
  imis_fit$lhood = lhood(imis_fit$resample)
  return(imis_fit)
}

plot_fit = function(imis_fit, pop_data, svy_data) {
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
  ggsave("temp.tiff", compression="lzw", dpi=600, width=2*3.42, height=2*2.44)
}

iso_code = "ZAF"
pop_data = readRDS("data/wpp-pop-male.rds")
pop_data = pop_data[pop_data$ISO_Alpha_3==iso_code,]
pop_data = pop_data[pop_data$Year >= 1970 & pop_data$Year <= 2025,]
svy_data = read.csv("data/survey-data.csv")
svy_data = svy_data[svy_data$ISO_Alpha_3==iso_code,]

imis_fit = fit_model(pop_data, svy_data)
plot_fit(imis_fit, pop_data, svy_data)


# ## +=+ BEGIN OLD CODE +=+
#
#
# X = log_imis(prior, likelihood, sample_prior, 500, 1000, 1000)
# # X = list(resample = sample.prior(25))
#
# post_par = list()
# for (k in 1:nrow(X$resample)) {
#   post_par[[k]] = list(mc_uptake_1 = X$resample[k, 1:6 ],
#                        mc_agedst_1 = X$resample[k, 7:8 ],
#                        mc_uptake_2 = X$resample[k, 9:12],
#                        mc_agedst_2 = X$resample[k,13:14])
# }
# imis_out = lapply(post_par, function(par) {model_sim(par, age_data)})
#
# X$prior = prior(X$resample)
# X$lhood = likelihood(X$resample)
# post_mode = which.max(X$prior + X$lhood)
#
# core = runif(nrow(X$resample), 0.7, 0.8)
# gset = rgb(core, core, core)
#
# plot.trend = function(a_lim, t_lim, col, main) {
#   mc_sset = subset(mc_data, age_min==a_lim[1] & age_max==a_lim[2] + 1)
#
#   a_ind = 1 + a_lim[1]:a_lim[2]
#   par(las=1, mar=c(2.5,3.0,1.0,1.0), cex=1)
#   plot(t_lim, c(0,1), cex=0, ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
#   abline(h=axTicks(2), lty=3, lwd=0.5)
#   for (k in 1:length(imis_out)) {with(imis_out[[k]], lines(year, rowSums(pop_crc[,a_ind]) / rowSums(pop_sum[,a_ind]), col=gset[k]))}
#   with(imis_out[[post_mode]], lines(year, rowSums(pop_crc[,a_ind]) / rowSums(pop_sum[,a_ind]), col=col))
#   points(mc_sset$year, mc_sset$prop_circumcised, pch=16, col=col)
#   arrows(mc_sset$year, mc_sset$ci_lower, mc_sset$year, mc_sset$ci_upper, code=3, angle=90, length=0.02, col=col)
#   axis(1, mgp=c(0.1, 0.3, 0.0), tck=-0.01, at=axTicks(1), labels=axTicks(1))
#   axis(2, mgp=c(0.1, 0.3, 0.0), tck=-0.01, at=axTicks(2), labels=100*axTicks(2))
#   title(ylab="Male circumcision prevalence, %", line=1.75)
#   title(main=main)
# }
#
# cset = brewer.pal(8,'Dark2')
# panel.names = sprintf("%d-%d", seq(15,50,5), seq(19,54,5))
# tiff(sprintf("%s-mc-data-2000-2030.tiff", tolower(isocode)), w=2*3.42, h=2*2.44, units="in", pointsize=8, compression="lzw", res=300)
# layout(matrix(1:8, nrow=2, byrow=TRUE))
# sapply(0:7, function(k) {plot.trend(c(15, 19) + 5 * k, c(2000, 2030), cset[k+1], panel.names[k+1])})
# dev.off()
#
# tiff(sprintf("%s-mc-post.tiff", tolower(isocode)), w=2*3.42, h=2*2.44, units="in", pointsize=8, compression="lzw", res=300)
# layout(matrix(1:4, nrow=2, byrow=TRUE))
# par(las=1, mar=c(2.5,3.0,1.0,1.0), cex=1)
# for (panel in 1:2) {
#   plot(c(1970,2020), c(0,1), cex=0, ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
#   for (i in 1:nrow(X$resample)) {lines(1970:2030, 1-exp(-mc_model_rate(1970:2030, post_par[[i]])[,panel]), col=gset[i])}
#   lines(1970:2030, 1-exp(-mc_model_rate(1970:2030, post_par[[post_mode]])[,panel]), col=cset[panel], lwd=1.5)
#   axis(1, mgp=c(0.1, 0.3, 0.0), tck=-0.01, at=axTicks(1), labels=axTicks(1))
#   axis(2, mgp=c(0.1, 0.3, 0.0), tck=-0.01, at=axTicks(2), labels=100*axTicks(2))
#   title(xlab="Year", line=1.5)
#   title(ylab=sprintf("Uptake trend %0.0f (%%/year)", panel), line=1.75)
# }
#
# for (panel in 1:2) {
#   ages = 0:54
#   plot(c(0,55), c(0,0.2), cex=0, ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
#   for (i in 1:nrow(X$resample)) {lines(ages, mc_model_dist(ages, post_par[[i]])[,panel], col=gset[i])}
#   lines(ages, mc_model_dist(ages, post_par[[post_mode]])[,panel], col=cset[panel], lwd=1.5)
#   axis(1, mgp=c(0.1, 0.3, 0.0), tck=-0.01, at=axTicks(1), labels=axTicks(1))
#   axis(2, mgp=c(0.1, 0.3, 0.0), tck=-0.01, at=axTicks(2), labels=sprintf("%0.1f", axTicks(2)))
#   title(xlab="Age", line=1.5)
#   title(ylab=sprintf("Uptake distribution %0.0f", panel), line=1.75)
# }
# dev.off()
#
# ## write out trends that can be input into the ASM
# post_out = model_sim(post_par[[post_mode]], age_data)
# post_cvg = with(post_out, cbind(rowSums(pop_crc[,15:19+1]) / rowSums(pop_sum[,15:19+1]),
#                                 rowSums(pop_crc[,20:24+1]) / rowSums(pop_sum[,20:24+1]),
#                                 rowSums(pop_crc[,25:29+1]) / rowSums(pop_sum[,25:29+1]),
#                                 rowSums(pop_crc[,30:34+1]) / rowSums(pop_sum[,30:34+1]),
#                                 rowSums(pop_crc[,35:39+1]) / rowSums(pop_sum[,35:39+1]),
#                                 rowSums(pop_crc[,40:44+1]) / rowSums(pop_sum[,40:44+1]),
#                                 rowSums(pop_crc[,45:49+1]) / rowSums(pop_sum[,45:49+1]),
#                                 rowSums(pop_crc[,50:54+1]) / rowSums(pop_sum[,50:54+1])))
# rownames(post_cvg) = imis_out[[post_mode]]$year
# colnames(post_cvg) = panel.names
# write.csv(t(post_cvg), sprintf('mc-trend-%s.csv', tolower(isocode)))\


