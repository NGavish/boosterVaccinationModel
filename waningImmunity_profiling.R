
# to be run after 'fitting' procedure


calcProfile <- function(par.ind,step,best.nll,best.est) {
  
  isProfileRun <<- TRUE
  profile <- list()
  best.par.est <- best.est[par.ind]
  par.name <- names(best.par.est)
  profile$par.vals <- best.par.est
  profile$nll.vals <- best.nll
  profile$par.est <- list()
  profile$par.est[[1]] <- best.est
  ii <- 2
  profile.par.ind <<- par.ind
  
  
  cur.nll <- best.nll
  par.est <<- best.est 
  profile.par.val <<- best.est[par.ind]
  while(cur.nll <= best.nll+2.75) {
    profile.par.val <<- profile.par.val + step
    source('waningImmunity_fitting.R')
    cur.nll <- best$value
    par.est <<- c(profile.par.val,best$par)[names(best.est)]
    profile$par.vals[ii] <- profile.par.val
    profile$nll.vals[ii] <- cur.nll
    profile$par.est[[ii]] <- par.est
    ii <- ii+1
  }
  
  cur.nll <- best.nll
  par.est <<- best.est
  profile.par.val <<- best.est[par.ind]
  while(cur.nll <= best.nll+2.75) {
    profile.par.val <<- profile.par.val - step
    source('waningImmunity_fitting.R')
    cur.nll <- best$value
    par.est <<- c(profile.par.val,best$par)[names(best.est)]
    profile$par.vals[ii] <- profile.par.val
    profile$nll.vals[ii] <- cur.nll
    profile$par.est[[ii]] <- par.est
    ii <- ii+1
  }
  ord.vals <- order(profile$par.vals)
  profile$par.vals <- profile$par.vals[ord.vals]
  profile$nll.vals <- profile$nll.vals[ord.vals]
  profile$par.est <- profile$par.est[ord.vals]
  names(profile$par.vals) <- rep(par.name,length(profile$par.vals))
  
  isProfileRun <<- FALSE
  return (profile)
}

fixProfile <- function(profile,par.ind,best.est) {
  isProfileRun <<- TRUE
  best.par.est <- best.est[par.ind]
  profile.par.ind <<- par.ind
  
  jj <- which(profile$par.vals==best.par.est)
  for(ii in jj:2) {
    profile.par.val <<- profile$par.vals[ii]
    nll <- profile$nll.vals[ii]
    if(nll > profile$nll.vals[ii-1]) {
      par.est <<- c(profile.par.val,profile$par.est[[ii-1]])[names(best.est)]
      source('waningImmunity_fitting.R')
      cur.nll <- best$value
      if(cur.nll < nll) {
        par.est <<- c(profile.par.val,best$par)[names(best.est)]
        profile$par.vals[ii] <- profile.par.val
        profile$nll.vals[ii] <- cur.nll
        profile$par.est[[ii]] <- par.est
      }
    }
  }
  for(ii in jj:(length(profile$par.vals)-1)) {
    profile.par.val <<- profile$par.vals[ii]
    nll <- profile$nll.vals[ii]
    if(nll > profile$nll.vals[ii+1]) {
      par.est <<- c(profile.par.val,profile$par.est[[ii+1]])[names(best.est)]
      source('waningImmunity_fitting.R')
      cur.nll <- best$value
      if(cur.nll < nll) {
        par.est <<- c(profile.par.val,best$par)[names(best.est)]
        profile$par.vals[ii] <- profile.par.val
        profile$nll.vals[ii] <- cur.nll
        profile$par.est[[ii]] <- par.est
      }
    }
  }
  isProfileRun <<- FALSE
  return (profile)
}

calcProfiles <- FALSE
if(calcProfiles) {
  par.num <- length(best.est)
  par.steps <- c(0.01,0.01,0.005,0.005,0.1,0.05,0.1,0.1)
  profiles <- list()
  for(par.ind in 1:par.num) {
    profiles[[par.ind]] <- calcProfile(par.ind,par.steps[par.ind],best.nll.val,best.est)
  }
  saveRDS(profiles,'likelihood_profiles.RDS')
  
  ci95  <- list()
  jj <- 0
  for(par.ind in 1:par.num) {
    best.par.est <- best.est[par.ind]
    ii <- which(profiles[[par.ind]]$nll.vals>best.nll.val+1.92 & profiles[[par.ind]]$par.vals<best.par.est)
    ci95[[jj+1]] <- profiles[[par.ind]]$par.est[[ii[length(ii)]]]
    ii <- which(profiles[[par.ind]]$nll.vals>best.nll.val+1.92 & profiles[[par.ind]]$par.vals>best.par.est)
    ci95[[jj+2]] <- profiles[[par.ind]]$par.est[[ii[1]]]
    jj <- jj+2 
  }
  sapply(1:length(ci95), function(i) round(ci95[[i]],2))
  saveRDS(ci95,'ci95_par_est.RDS')
}


ci95 <- readRDS('ci95_par_est.RDS')
profiles <- readRDS('likelihood_profiles.RDS')
plotProfiles <- TRUE
if(plotProfiles) {
  par.num <- length(best.est)
  par.plot.names <- c(expression(R[t[0]]),expression(delta[1]),
                      expression(rho[1]),expression(rho[9]),
                      expression(omega[w]),expression(omega[s]),expression(omega[c]),'r')
  x11(width=8,height=10); par(mfrow=c(4,2),mar=c(4,5,2,4))
  for(par.ind in 1:par.num) {
    curve1 <- data.frame(x=profiles[[par.ind]]$par.vals,y=profiles[[par.ind]]$nll.vals)
    curve2 <- data.frame(x=profiles[[par.ind]]$par.vals,y=best.nll.val+1.92)
    plot(curve1$x, curve1$y,xlab=par.plot.names[par.ind],ylab='Negative log-likelihood',
         cex.lab=1.25,cex.axis=1.25,ylim=c(15136,15140))
    lines(curve1$x, curve1$y,lwd=2)
    abline(h=curve2$y,col='red',lwd=2)
  }
}


