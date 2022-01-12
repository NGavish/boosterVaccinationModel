
fitModel <- function(par.est, parms, 
                     fit_data, fit_period, 
                     par.min=NA, par.max=NA, ...) {
  
  if (missing(par.est)) {
    stop('missing par.est')
  }
  
  par.est.names <- names(par.est)
  if(!is.na(par.min[1]))
    par.est <- sapply(1:length(par.est), function(i) max(par.est[i],par.min[i]))
  if(!is.na(par.max[1]))
    par.est <- sapply(1:length(par.est), function(i) min(par.est[i],par.max[i]))
  names(par.est) <- par.est.names
  
  parms <- extractParameters(par.est, parms)
  Cij <- computeCij(parms)
  
  l1 <- remove_hidden_infections_from_schedules(parms,vaccinationSchedule0,boostSchedule0)
  vaccinationSchedule <- l1$vaccinationSchedule
  boostSchedule <- l1$boostSchedule
  
  l2 <- set_initial_conditions(parms,vaccinationSchedule)
  i0_nv <- l2$i0_nv
  i0_v  <- l2$i0_v
  S0_nv <- l2$S0_nv
  S0_v  <- l2$S0_v
  
  parms$Tmax <- length(fit_period)
  model_out <- run_ag_aoi_waningImmunity(parms,S0_v,S0_nv,i0_nv,i0_v,Cij,vaccinationSchedule,boostSchedule)
  
  loss.vals <- rep(0,6)
  if(fitDetected) {
    loss.vals[1] <- sum(sapply(1:m,function(i) dpois(fit_data$detected_nv[fit_period,i],lambda=(model_out$detected_nv[fit_period,i]+1e-5),log=T)))
    loss.vals[2] <- sum(sapply(1:m,function(i) dpois(fit_data$detected_v[fit_period,i],lambda=(model_out$detected_v[fit_period,i]+1e-5),log=T)))
    loss.vals[3] <- sum(sapply(1:m,function(i) dpois(fit_data$detected_b[fit_period,i],lambda=(model_out$detected_b[fit_period,i]+1e-5),log=T)))
  }
  if(fitSevere) {
    
    loss.vals[4] <- sum(sapply(1:m,function(i) dpois(fit_data$severe_nv[fit_period,i],lambda=(model_out$sevp_nv[fit_period,i]+1e-5),log=T)))
    loss.vals[5] <- sum(sapply(1:m,function(i) dpois(fit_data$severe_v[fit_period,i],lambda=(model_out$sevp_v[fit_period,i]+1e-5),log=T)))
    loss.vals[6] <- sum(sapply(1:m,function(i) dpois(fit_data$severe_b[fit_period,i],lambda=(model_out$sevp_b[fit_period,i]+1e-5),log=T)))
  }
  loss <- -sum(loss.vals)
  
  # if(is.nan(loss)) {
  #   print('error')
  # }
  # print(par.est)
  print(loss)
  
  return (loss)
}


extractParameters <- function(par.est, parms) {
  
  if(isProfileRun) {
    par.est <- c(par.est,profile.par.val)
  }

  parms$beta_multiplier <- par.est[beta.mult.names]

  parms$susceptibilityFactor[1] <- par.est[child.sus.factor.names]
  
  dr <- par.est[det.rates.names]
  parms$det_rates_nv <- c(rep(dr[1],7),rep(dr[2],2))
  parms$det_rates_v <-  parms$det_rates_nv
  
  return (parms)
}

computeCij <- function(parms) {
  Cij <- lapply(1:parms$Tmax, function(t)
                mobility_data$household[t]*Cij_household+
                mobility_data$work[t]*Cij_work+
                mobility_data$school[t]*Cij_school+
                mobility_data$community[t]*Cij_community)
  
  Cij <- lapply(1:length(Cij), function(t) Cij[[t]]*parms$susceptibilityFactor)
  Cij <- lapply(1:length(Cij), function(t) t(t(Cij[[t]])*parms$infectivityFactor))
  norm_fact <- parms$beta_multiplier/max(Re(eigen(Cij[[1]])$values))
  Cij <- lapply(1:length(Cij), function(t) Cij[[t]]*norm_fact)
  return(Cij)
}

##################################################
# Fitting procedure
##################################################

if(!exists('isProfileRun'))
  isProfileRun <- FALSE

beta.mult.names <- 'R'
child.sus.factor.names <- 'delta1'
det.rates.names <- paste0('rho',1:2)

#R
par.guess <- 2.5
par.min <- 1
par.max <- 10
par.names <- beta.mult.names
               
#sus factor
par.guess <- c(par.guess,rep(0.5,length(child.sus.factor.names)))
par.min <- c(par.min,rep(0,length(child.sus.factor.names)))
par.max <- c(par.max,rep(1,length(child.sus.factor.names)))
par.names <- c(par.names,child.sus.factor.names) 

#det rates
par.guess <- c(par.guess,rep(0.5,length(det.rates.names)))
par.min <- c(par.min,rep(0.25,length(det.rates.names)))
par.max <- c(par.max,rep(1,length(det.rates.names)))
par.names <- c(par.names,det.rates.names)

par.est <- c(2.605,0.601,0.446,0.572) #loss=44450
# par.est <- par.guess
names(par.est) <- par.names

doFit <- FALSE
fitDetected <- TRUE
fitSevere <- FALSE

if(isProfileRun) {
  par.est <- par.est[-profile.par.ind]
  par.min <- par.min[-profile.par.ind]
  par.max <- par.max[-profile.par.ind]
  par.names <- par.names[-profile.par.ind]
} 

if(doFit) {
  for(jj in 1:3) {
    best <- optim(par.est, fitModel, method='Nelder-Mead', #'BFGS', #
                  control=list(parscale=pmax(1,par.est/min(par.est[which(par.est>0)])),maxit=200,reltol=1e-4),
                  parms=parms, fit_data=fit_data, fit_period=fit_period, par.min=par.min, par.max=par.max)
    par.est <- best$par
  }
  print(paste0('loss: ',round(best$value)))
}

par.est <- sapply(1:length(par.est), function(i) max(par.est[i],par.min[i]))
par.est <- sapply(1:length(par.est), function(i) min(par.est[i],par.max[i]))
names(par.est) <- par.names
parms <- extractParameters(par.est, parms)
Cij <- computeCij(parms)

l1 <- remove_hidden_infections_from_schedules(parms,vaccinationSchedule0,boostSchedule0)
vaccinationSchedule <- l1$vaccinationSchedule
boostSchedule <- l1$boostSchedule

l2 <- set_initial_conditions(parms,vaccinationSchedule)
i0_nv <- l2$i0_nv
i0_v  <- l2$i0_v
S0_nv <- l2$S0_nv
S0_v  <- l2$S0_v

print(round(par.est,3))

model_out <-  run_ag_aoi_waningImmunity(parms,S0_v,S0_nv,i0_nv,i0_v,Cij,vaccinationSchedule,boostSchedule)




