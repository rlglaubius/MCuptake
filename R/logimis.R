## Update of Raftery and Bao's IMIS R implementation to use log-transformed priors and likelihoods
log_imis = function (prior, likelihood, sample.prior, B = 1000, B.re = 3000, number_k = 100, D = 0) {
  B0 = B * 10
  X_all = X_k = sample.prior(B0)
  if (is.vector(X_all))
    Sig2_global = var(X_all)
  if (is.matrix(X_all))
    Sig2_global = cov(X_all)
  stat_all = matrix(NA, 6, number_k)
  center_all = prior_all = like_all = NULL
  sigma_all = list()
  if (D >= 1)
    option.opt = 1
  if (D == 0) {
    option.opt = 0
    D = 1
  }
  for (k in 1:number_k) {
    ptm.like = proc.time()
    prior_all = c(prior_all, prior(X_k))
    like_all = c(like_all, likelihood(X_k))
    ptm.use = (proc.time() - ptm.like)[3]
    if (k == 1)
      print(paste(B0, "likelihoods are evaluated in", round(ptm.use/60,
                                                            2), "minutes"))
    if (k == 1)
      envelop_all = prior_all
    if (k > 1)
      envelop_all = log(apply(rbind(exp(prior_all) * B0/B, gaussian_all),
                          2, sum)/(B0/B + D + (k - 2)))
    Weights = exp(prior_all + like_all - envelop_all - max(like_all))
    stat_all[1, k] = log(mean(exp(prior_all + like_all - envelop_all)))
    Weights = Weights/sum(Weights)
    stat_all[2, k] = sum(1 - (1 - Weights)^B.re)
    stat_all[3, k] = max(Weights)
    stat_all[4, k] = 1/sum(Weights^2)
    stat_all[5, k] = -sum(Weights * log(Weights), na.rm = TRUE)/log(length(Weights))
    stat_all[6, k] = var(Weights/mean(Weights))
    if (k == 1)
      print("Stage   MargLike   UniquePoint   MaxWeight   ESS")
    print(c(k, round(stat_all[1:4, k], 3)))
    if (k == 1 & option.opt == 1) {
      if (is.matrix(X_all))
        Sig2_global = cov(X_all[which(like_all > min(like_all)),
                                ])
      X_k = which_exclude = NULL
      label_weight = sort(Weights, decreasing = TRUE, index = TRUE)
      which_remain = which(Weights > label_weight$x[B0])
      size_remain = length(which_remain)
      for (i in 1:D) {
        important = NULL
        if (length(which_remain) > 0)
          important = which_remain[which(Weights[which_remain] ==
                                           max(Weights[which_remain]))]
        if (length(important) > 1)
          important = sample(important, 1)
        if (is.vector(X_all))
          X_imp = X_all[important]
        if (is.matrix(X_all))
          X_imp = X_all[important, ]
        which_exclude = union(which_exclude, important)
        which_remain = setdiff(which_remain, which_exclude)
        posterior = function(theta) {
          -prior(theta) - likelihood(theta)
        }
        if (is.vector(X_all)) {
          if (length(important) == 0)
            X_imp = center_all[1]
          optimizer = optim(X_imp, posterior, method = "BFGS",
                            hessian = TRUE, control = list(parscale = sqrt(Sig2_global)/10,
                                                           maxit = 5000))
          print(paste("maximum posterior=", round(-optimizer$value,
                                                  2), ", likelihood=", round(likelihood(optimizer$par),
                                                                             2), ", prior=", round(prior(optimizer$par),
                                                                                                   2), ", time used=", round(ptm.use/60, 2),
                      "minutes, convergence=", optimizer$convergence))
          center_all = c(center_all, optimizer$par)
          sigma_all[[i]] = solve(optimizer$hessian)
          X_k = c(X_k, rnorm(B, optimizer$par, sqrt(sigma_all[[i]])))
          distance_remain = abs(X_all[which_remain] -
                                  optimizer$par)
        }
        if (is.matrix(X_all)) {
          if (length(important) == 0)
            X_imp = center_all[1, ]
          ptm.opt = proc.time()
          optimizer = optim(X_imp, posterior, method = "Nelder-Mead",
                            control = list(maxit = 10000, parscale = sqrt(diag(Sig2_global))))
          theta.NM = optimizer$par
          optimizer = optim(theta.NM, posterior, method = "BFGS",
                            hessian = TRUE, control = list(parscale = sqrt(diag(Sig2_global)),
                                                           maxit = 1000, ndeps=rep(1e-6, length(theta.NM))))
          ptm.use = (proc.time() - ptm.opt)[3]
          print(paste("maximum posterior=", round(-optimizer$value,
                                                  2), ", likelihood=", round(likelihood(optimizer$par),
                                                                             2), ", prior=", round(prior(optimizer$par),
                                                                                                   2), ", time used=", round(ptm.use/60, 2),
                      "minutes, convergence=", optimizer$convergence))
          center_all = rbind(center_all, optimizer$par)
          if (min(eigen(optimizer$hessian)$values) >
              0)
            sigma_all[[i]] = solve(optimizer$hessian)
          if (min(eigen(optimizer$hessian)$values) <=
              0) {
            eigen.values = eigen(optimizer$hessian)$values
            eigen.values[which(eigen.values < 0)] = 0
            hessian = eigen(optimizer$hessian)$vectors %*%
              diag(eigen.values) %*% t(eigen(optimizer$hessian)$vectors)
            sigma_all[[i]] = solve(hessian + diag(1/diag(Sig2_global)))
          }
          X_k = rbind(X_k, mvtnorm::rmvnorm(B, optimizer$par,
                                   sigma_all[[i]]))
          # distance_remain = mahalanobis(X_all[which_remain,
          #                                     ], optimizer$par, diag(diag(Sig2_global)))
          distance_remain = mahalanobis(X_all[which_remain,], optimizer$par, Sig2_global)
        }
        label_dist = sort(distance_remain, decreasing = FALSE,
                          index = TRUE)
        which_exclude = union(which_exclude, which_remain[label_dist$ix[1:floor(size_remain/D)]])
        which_remain = setdiff(which_remain, which_exclude)
      }
      if (is.matrix(X_all))
        X_all = rbind(X_all, X_k)
      if (is.vector(X_all))
        X_all = c(X_all, X_k)
    }
    if (k > 1 | option.opt == 0) {
      important = which(Weights == max(Weights))
      if (length(important) > 1)
        important = important[1]
      if (is.matrix(X_all))
        X_imp = X_all[important, ]
      if (is.vector(X_all))
        X_imp = X_all[important]
      if (is.matrix(X_all))
        center_all = rbind(center_all, X_imp)
      if (is.vector(X_all))
        center_all = c(center_all, X_imp)
      if (is.matrix(X_all))
        # distance_all = mahalanobis(X_all, X_imp, diag(diag(Sig2_global)))
        distance_all = mahalanobis(X_all, X_imp, Sig2_global)
      if (is.vector(X_all))
        distance_all = abs(X_all - X_imp)
      label_nr = sort(distance_all, decreasing = FALSE,
                      index = TRUE)
      which_var = label_nr$ix[1:B]
      if (is.matrix(X_all))
        Sig2 = cov.wt(X_all[which_var, ], wt = Weights[which_var] +
                        1/length(Weights), cor = FALSE, center = X_imp,
                      method = "unbias")$cov
      if (is.vector(X_all)) {
        Weights_var = Weights[which_var] + 1/length(X_all)
        Weights_var = Weights_var/sum(Weights_var)
        Sig2 = (X_all[which_var] - X_imp)^2 %*% Weights_var
      }
      sigma_all[[D + k - 1]] = Sig2
      if (is.matrix(X_all))
        X_k = mvtnorm::rmvnorm(B, X_imp, Sig2)
      if (is.vector(X_all))
        X_k = rnorm(B, X_imp, sqrt(Sig2))
      if (is.matrix(X_all))
        X_all = rbind(X_all, X_k)
      if (is.vector(X_all))
        X_all = c(X_all, X_k)
    }
    if (k == 1) {
      gaussian_all = matrix(NA, D, B0 + D * B)
      for (i in 1:D) {
        if (is.matrix(X_all))
          gaussian_all[i, ] = mvtnorm::dmvnorm(X_all, center_all[i,
                                                        ], sigma_all[[i]])
        if (is.vector(X_all))
          gaussian_all[i, ] = dnorm(X_all, center_all[i],
                                    sqrt(sigma_all[[i]]))
      }
    }
    if (k > 1) {
      if (is.vector(X_all))
        gaussian_new = matrix(0, D + k - 1, length(X_all))
      if (is.matrix(X_all))
        gaussian_new = matrix(0, D + k - 1, dim(X_all)[1])
      if (is.matrix(X_all)) {
        gaussian_new[1:(D + k - 2), 1:(dim(X_all)[1] -
                                         B)] = gaussian_all
        gaussian_new[D + k - 1, ] = mvtnorm::dmvnorm(X_all, X_imp,
                                            sigma_all[[D + k - 1]])
        for (j in 1:(D + k - 2)) gaussian_new[j, (dim(X_all)[1] -
                                                    B + 1):dim(X_all)[1]] = mvtnorm::dmvnorm(X_k, center_all[j,
                                                                                                    ], sigma_all[[j]])
      }
      if (is.vector(X_all)) {
        gaussian_new[1:(D + k - 2), 1:(length(X_all) -
                                         B)] = gaussian_all
        gaussian_new[D + k - 1, ] = dnorm(X_all, X_imp,
                                          sqrt(sigma_all[[D + k - 1]]))
        for (j in 1:(D + k - 2)) gaussian_new[j, (length(X_all) -
                                                    B + 1):length(X_all)] = dnorm(X_k, center_all[j],
                                                                                  sqrt(sigma_all[[j]]))
      }
      gaussian_all = gaussian_new
    }
    if (stat_all[2, k] > (1 - exp(-1)) * B.re)
      break
  }
  nonzero = which(Weights > 0)
  which_X = sample(nonzero, B.re, replace = TRUE, prob = Weights[nonzero])
  if (is.matrix(X_all))
    resample_X = X_all[which_X, ]
  if (is.vector(X_all))
    resample_X = X_all[which_X]
  return(list(stat = t(stat_all), resample = resample_X, center = center_all))
}
