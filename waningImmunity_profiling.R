
# to be run after 'fitting' procedure

mle_best <- best
best.nll.val <- mle_best$value
best.par.est <- mle_best$par

calcProfiles <- FALSE

if(calcProfiles) {
  
  isProfileRun <- TRUE

  profiles <- list()
  
  profile.par.vals <- list()
  profile.par.vals[[1]] <- seq(2.600,2.610,0.002)
  profile.par.vals[[2]] <- seq(0.595,0.605,0.002)
  profile.par.vals[[3]] <- seq(0.443,0.450,0.001)
  profile.par.vals[[4]] <- seq(0.560,0.590,0.005)

  for(profile.par.ind in 1:length(profile.par.vals)) {
    
    profile.nll.vals <- c()
    profile.est.vals <- list()
    ii <- 1
    for(profile.par.val in profile.par.vals[[profile.par.ind]]) {
      profile.par.name <- names(best.par.est[profile.par.ind])
      names(profile.par.val) <- profile.par.name
      source('waningImmunity_fitting.R')
      profile.nll.vals[ii] <- best$value
      profile.est.vals[[ii]] <- par.est
      ii<- ii+1
    }
    profiles[[profile.par.ind]] <- list(par.vals=profile.par.vals[[profile.par.ind]],nll.vals=profile.nll.vals,est.vals=profile.est.vals)
  }
  
  profiles[[1]]$ci95 <- c(2.6028,2.6072)
  profiles[[2]]$ci95 <- c(0.5982,0.6036)
  profiles[[3]]$ci95 <- c(0.4452,0.4480)
  profiles[[4]]$ci95 <- c(0.5655,0.5821)
  for(profile.par.ind in 1:length(profiles)) {
    profiles[[profile.par.ind]]$ci95_model_out  <- list()
    profiles[[profile.par.ind]]$ci95_par_est  <- list()
    profile.par.name <- names(best.par.est[profile.par.ind])
    for(ii in 1:2) {
      profile.par.val <- profiles[[profile.par.ind]]$ci95[ii]
      names(profile.par.val) <- profile.par.name
      source('waningImmunity_fitting.R')
      profiles[[profile.par.ind]]$ci95_model_out[[ii]] <- model_out
      profiles[[profile.par.ind]]$ci95_par_est[[ii]] <- c(par.est,profile.par.val)[names(best.par.est)]
    }
  }
  names(profiles) <- names(best.par.est)
  isProfileRun <- FALSE
  best <- mle_best
  saveRDS(profiles,'likelihood_profiles.RDS')
  
  ci95_model_out <- list()
  ci95_par_est <- list()
  for(i in 1:length(profiles)) {
    ci95_model_out[[(i-1)*2+1]] <- profiles[[i]]$ci95_model_out[[1]]
    ci95_model_out[[(i-1)*2+2]] <- profiles[[i]]$ci95_model_out[[2]]
    ci95_par_est[[(i-1)*2+1]] <- profiles[[i]]$ci95_par_est[[1]]
    ci95_par_est[[(i-1)*2+2]] <- profiles[[i]]$ci95_par_est[[2]]
  }
  saveRDS(ci95_par_est,'ci95_par_est.RDS')
}

profiles <- readRDS('likelihood_profiles.RDS')
plotProfiles <- TRUE
if(plotProfiles) {
  par.plot.names <- c(expression('R(t_0)'),expression(delta~1),
                      expression(rho~ '1:7'),expression(rho~ '8:9'))
  x11(width=8,height=8); par(mfrow=c(2,2))
  for(profile.par.ind in 1:length(profiles)) {
    curve1 <- data.frame(x=profiles[[profile.par.ind]]$par.vals,y=profiles[[profile.par.ind]]$nll.vals)
    curve2 <- data.frame(x=profiles[[profile.par.ind]]$par.vals,y=best.nll.val+1.92)
    plot(curve1$x, curve1$y,xlab=par.plot.names[profile.par.ind],ylab='Negative Log-Likelihood',ylim=c(44448,44460),size=3,cex.lab=1.5,cex.axis=1.5)
    lines(curve1$x, curve1$y,lwd=2)
    abline(h=curve2$y,col='red',lwd=2)
    abline(v=profiles[[profile.par.ind]]$ci95[1],lty=2,lwd=2)
    abline(v=profiles[[profile.par.ind]]$ci95[2],lty=2,lwd=2)
  }
}



