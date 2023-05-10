library(RColorBrewer)
library(signal)
source("R/logimis.R")

isocode = "LSO"

age_data = readRDS("data/wpp-pop-male.rds")
age_data = age_data[age_data$ISO_Alpha_3==isocode & age_data$Year>=1970,]

mc_data = read.csv(sprintf("data/mc-data-%s.csv", tolower(isocode)))

par = list(mc_uptake_1 = c(-1.0, 1.0, 0.0, 0.0, 1.0, 1.0),
           mc_uptake_2 = c( 0.0, 0.0, 0.0, 0.0),
           mc_agedst_1 = c(3.0, 1.0),
           mc_agedst_2 = c(3.0, 1.0))

## Male circumcision uptake rates over time
mc_model_rate = function(year, par) {
  z = par$mc_uptake_1[1] + (par$mc_uptake_1[2] - par$mc_uptake_1[1]) / (1.0 + exp(-par$mc_uptake_1[3] * (year - 1970 - par$mc_uptake_1[4])))
  mc_rate_1 = exp(-0.5 * (z / par$mc_uptake_1[5])^2) * par$mc_uptake_1[6]
  mc_rate_2 = par$mc_uptake_2[1] + (par$mc_uptake_2[2] - par$mc_uptake_2[1]) / (1.0 + exp(-par$mc_uptake_2[3] * (year - 1970 - par$mc_uptake_2[4])))
  return(cbind(mc_rate_1, mc_rate_2))
}

## Age distribution of male circumcision
mc_model_dist = function(age, par) {
  return(cbind(dnbinom(age, size=par$mc_agedst_1[1], mu=par$mc_agedst_1[2]),
               dnbinom(age, size=par$mc_agedst_2[1], mu=par$mc_agedst_2[2])))
}

mc_model = function(year, age, par) {
  mc_rate = mc_model_rate(year, par)
  mc_dist = mc_model_rate(age,  par)
  return(outer(mc_rate[,1], mc_dist[,1]) + outer(mc_rate[,2], mc_dist[,2]))
}

model_sim = function(par, age_data) {
  year = unique(age_data$Year)
  ages = unique(age_data$Age)
  num_a = length(ages)
  num_t = length(year)

  age_wide = reshape2::dcast(age_data, Year~Age, value.var="Value")
  pop_sum = data.matrix(age_wide[,2:ncol(age_wide)])
  pop_unc = matrix(0, ncol=num_a, nrow=num_t) # uncircumcised population
  pop_crc = matrix(0, ncol=num_a, nrow=num_t) # circumcised population

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

  ## Initialize prevalence at age 0 ()
  pop_crc[,1] = pop_sum[,1] * mc_prob[,1]
  pop_unc[,1] = pop_sum[,1] - pop_crc[,1]

  ## Calculate MC prevalence from uptake and age 0 and base-year levels
  for (yi in 2:num_t) {
    sx = pop_sum[yi,2:num_a] / pop_sum[yi-1,1:(num_a-1)] # approximate survival from year-to-year populations
    pop_crc[yi,2:num_a] = (pop_crc[yi-1,1:(num_a-1)] + mc_prob[yi,2:num_a] * pop_unc[yi-1,1:(num_a-1)]) * sx
    pop_unc[yi,2:num_a] = pop_sum[yi,2:num_a] - pop_crc[yi,2:num_a]
  }

  return(list(year = year, pop_sum = pop_sum, pop_unc = pop_unc, pop_crc = pop_crc))
}

sample_prior = function(n) {
  X = cbind(rnorm(n, 0, 1),  ## mc rate, limit as t -> -inf in probit space
            rnorm(n, 0, 1),  ## mc rate, limit as t -> +inf in probit space
            rexp(n, 0.5),    ## mc rate, maximum slope magnitude (prior: less rapid changes are more likely)
            rnorm(n, 2012-1970, 8), ## mc rate, time of change
            rexp(n, 1.0),    ## mc rate maximum
            rexp(n, 1.0),    ## mc rate probit variance

            rexp(n, 0.1),    ## mc age distribution, negative binomial size
            runif(n, 0, 50), ## mc age distribution, negative binomial mean

            rexp(n, 8.0),
            rexp(n, 8.0),
            rexp(n, 0.5),
            rnorm(n, 2012-1970, 8),

            rexp(n, 0.1),
            runif(n, 0, 50))
}

prior = function(X) {
  return(rowSums(
    cbind(
      dnorm(X[,1], 0, 1,  log=TRUE),
      dnorm(X[,2], 0, 1,  log=TRUE),
      dexp( X[,3], 0.5,   log=TRUE),
      dnorm(X[,4], 2012-1970, 8, log=TRUE),
      dexp( X[,5], 1.0,   log=TRUE),
      dexp( X[,6], 1.0,   log=TRUE),

      dexp( X[,7], 0.1,  log=TRUE),
      dunif(X[,8], 0, 50,log=TRUE),

      dexp( X[, 9], 8.0,   log=TRUE),
      dexp( X[,10], 8.0,   log=TRUE),
      dexp( X[,11], 0.5,   log=TRUE),
      dnorm(X[,12], 2012-1970, 8, log=TRUE),

      dexp( X[,13], 0.1,  log=TRUE),
      dunif(X[,14], 0, 50,log=TRUE))))
}

likelihood = function(X) {
  P = prior(X)
  L = rep(-Inf, nrow(X))

  for (k in which(is.finite(P))) {
    inputs = par
    inputs$mc_uptake_1 = X[k, 1:6 ]
    inputs$mc_agedst_1 = X[k, 7:8 ]
    inputs$mc_uptake_2 = X[k, 9:12]
    inputs$mc_agedst_2 = X[k,13:14]
    output = model_sim(inputs, age_data)

    L[k] = sum(apply(mc_data[,c("age_min", "age_max", "year", "num_uncircumcised", "num_circumcised")], 1, function(obs) {
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

# X = IMIS.log(prior, likelihood, sample_prior, 500, 1000, 1000)
# X = list(resample = sample.prior(25))
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

plot.trend = function(a_lim, t_lim, col, main) {
  mc_sset = subset(mc_data, age_min==a_lim[1] & age_max==a_lim[2] + 1)

  a_ind = 1 + a_lim[1]:a_lim[2]
  par(las=1, mar=c(2.5,3.0,1.0,1.0), cex=1)
  plot(t_lim, c(0,1), cex=0, ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
  abline(h=axTicks(2), lty=3, lwd=0.5)
  for (k in 1:length(imis_out)) {with(imis_out[[k]], lines(year, rowSums(pop_crc[,a_ind]) / rowSums(pop_sum[,a_ind]), col=gset[k]))}
  with(imis_out[[post_mode]], lines(year, rowSums(pop_crc[,a_ind]) / rowSums(pop_sum[,a_ind]), col=col))
  points(mc_sset$year, mc_sset$prop_circumcised, pch=16, col=col)
  arrows(mc_sset$year, mc_sset$ci_lower, mc_sset$year, mc_sset$ci_upper, code=3, angle=90, length=0.02, col=col)
  axis(1, mgp=c(0.1, 0.3, 0.0), tck=-0.01, at=axTicks(1), labels=axTicks(1))
  axis(2, mgp=c(0.1, 0.3, 0.0), tck=-0.01, at=axTicks(2), labels=100*axTicks(2))
  title(ylab="Male circumcision prevalence, %", line=1.75)
  title(main=main)
}

cset = brewer.pal(8,'Dark2')
panel.names = sprintf("%d-%d", seq(15,50,5), seq(19,54,5))
tiff(sprintf("%s-mc-data-2000-2030.tiff", tolower(isocode)), w=2*3.42, h=2*2.44, units="in", pointsize=8, compression="lzw", res=300)
layout(matrix(1:8, nrow=2, byrow=TRUE))
sapply(0:7, function(k) {plot.trend(c(15, 19) + 5 * k, c(2000, 2030), cset[k+1], panel.names[k+1])})
dev.off()

# tiff(sprintf("%s-mc-data-2000-2020.tiff", tolower(country)), w=2*3.42, h=2*2.44, units="in", pointsize=8, compression="lzw", res=300)
# layout(matrix(1:8, nrow=2, byrow=TRUE))
# sapply(0:7, function(k) {plot.trend(c(15, 19) + 5 * k, c(2000, 2020), cset[k+1], panel.names[k+1])})
# dev.off()
#
# tiff(sprintf("%s-mc-post.tiff", tolower(country)), w=2*3.42, h=2*2.44, units="in", pointsize=8, compression="lzw", res=300)
# layout(matrix(1:4, nrow=2, byrow=TRUE))
# par(las=1, mar=c(2.5,3.0,1.0,1.0), cex=1)
# for (panel in 1:2) {
#   plot(c(1970,2020), c(0,1), cex=0, ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
#   for (i in 1:nrow(X$resample)) {lines(1970:2020, 1-exp(-mc.model.rate(1970:2020, post.par[[i]])[,panel]), col=gset[i])}
#   lines(1970:2020, 1-exp(-mc.model.rate(1970:2020, post.par[[post.mode]])[,panel]), col=cset[panel], lwd=1.5)
#   axis(1, mgp=c(0.1, 0.3, 0.0), tck=-0.01, at=axTicks(1), labels=axTicks(1))
#   axis(2, mgp=c(0.1, 0.3, 0.0), tck=-0.01, at=axTicks(2), labels=100*axTicks(2))
#   title(xlab="Year", line=1.5)
#   title(ylab=sprintf("Uptake trend %0.0f (%%/year)", panel), line=1.75)
# }
#
# for (panel in 1:2) {
#   ages = 0:54
#   plot(c(0,55), c(0,0.2), cex=0, ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
#   for (i in 1:nrow(X$resample)) {lines(ages, mc.model.dist(ages, post.par[[i]])[,panel], col=gset[i])}
#   lines(ages, mc.model.dist(ages, post.par[[post.mode]])[,panel], col=cset[panel], lwd=1.5)
#   axis(1, mgp=c(0.1, 0.3, 0.0), tck=-0.01, at=axTicks(1), labels=axTicks(1))
#   axis(2, mgp=c(0.1, 0.3, 0.0), tck=-0.01, at=axTicks(2), labels=sprintf("%0.1f", axTicks(2)))
#   title(xlab="Age", line=1.5)
#   title(ylab=sprintf("Uptake distribution %0.0f", panel), line=1.75)
# }
# dev.off()
#
# ## write out trends that can be input into the ASM
# post.out = model.sim(post.par[[post.mode]])
# post.cvg = with(post.out, cbind(rowSums(pop.crc[,15:19+1]) / rowSums(pop.sum[,15:19+1]),
#                                 rowSums(pop.crc[,20:24+1]) / rowSums(pop.sum[,20:24+1]),
#                                 rowSums(pop.crc[,25:29+1]) / rowSums(pop.sum[,25:29+1]),
#                                 rowSums(pop.crc[,30:34+1]) / rowSums(pop.sum[,30:34+1]),
#                                 rowSums(pop.crc[,35:39+1]) / rowSums(pop.sum[,35:39+1]),
#                                 rowSums(pop.crc[,40:44+1]) / rowSums(pop.sum[,40:44+1]),
#                                 rowSums(pop.crc[,45:49+1]) / rowSums(pop.sum[,45:49+1]),
#                                 rowSums(pop.crc[,50:54+1]) / rowSums(pop.sum[,50:54+1])))
# rownames(post.cvg) = age.data$year
# colnames(post.cvg) = panel.names
# write.csv(t(post.cvg), sprintf('mc-trend-%s.csv', tolower(country)))


