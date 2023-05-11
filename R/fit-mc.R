source("R/model-defs.R") # TODO: Remove
source("R/logimis.R")    # TODO: Remove

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
  imis_fit = log_imis(prior, lhood, sample_prior, 500, 1000, 1000)
  imis_fit$prior = prior(imis_fit$resample)
  imis_fit$lhood = lhood(imis_fit$resample)
  return(imis_fit)
}


# iso_code = "UGA"
# pop_data = readRDS("data/wpp-pop-male.rds")
# pop_data = pop_data[pop_data$ISO_Alpha_3==iso_code,]
# pop_data = pop_data[pop_data$Year >= 1970 & pop_data$Year <= 2025,]
# svy_data = read.csv("data/survey-data.csv")
# svy_data = svy_data[svy_data$ISO_Alpha_3==iso_code,]
#
# imis_fit = fit_model(pop_data, svy_data)
# plot_fit(imis_fit, pop_data, svy_data)


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


