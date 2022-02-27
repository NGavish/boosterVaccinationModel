
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
    if(fitUsingNB) {
      loss.vals[1] <- sum(sapply(1:m,function(i) dnbinom(fit_data$detected_nv[fit_period,i],mu=(model_out$detected_nv[fit_period,i]+1e-5),size=parms$dp,log=T)))
      loss.vals[2] <- sum(sapply(1:m,function(i) dnbinom(fit_data$detected_v[fit_period,i],mu=(model_out$detected_v[fit_period,i]+1e-5),size=parms$dp,log=T)))
      loss.vals[3] <- sum(sapply(1:m,function(i) dnbinom(fit_data$detected_b[fit_period,i],mu=(model_out$detected_b[fit_period,i]+1e-5),size=parms$dp,log=T)))
    } else {
      loss.vals[1] <- sum(sapply(1:m,function(i) dpois(fit_data$detected_nv[fit_period,i],lambda=(model_out$detected_nv[fit_period,i]+1e-5),log=T)))
      loss.vals[2] <- sum(sapply(1:m,function(i) dpois(fit_data$detected_v[fit_period,i],lambda=(model_out$detected_v[fit_period,i]+1e-5),log=T)))
      loss.vals[3] <- sum(sapply(1:m,function(i) dpois(fit_data$detected_b[fit_period,i],lambda=(model_out$detected_b[fit_period,i]+1e-5),log=T)))
    }
  }
  if(fitSevere) {
    if(fitUsingNB) {
      loss.vals[4] <- sum(sapply(1:m,function(i) dnbinom(fit_data$severe_nv[fit_period,i],mu=(model_out$sevp_nv[fit_period,i]+1e-5),size=parms$dp,log=T)))
      loss.vals[5] <- sum(sapply(1:m,function(i) dnbinom(fit_data$severe_v[fit_period,i],mu=(model_out$sevp_v[fit_period,i]+1e-5),size=parms$dp,log=T)))
      loss.vals[6] <- sum(sapply(1:m,function(i) dnbinom(fit_data$severe_b[fit_period,i],mu=(model_out$sevp_b[fit_period,i]+1e-5),size=parms$dp,log=T)))
    } else {
      loss.vals[4] <- sum(sapply(1:m,function(i) dpois(fit_data$severe_nv[fit_period,i],lambda=(model_out$sevp_nv[fit_period,i]+1e-5),log=T)))
      loss.vals[5] <- sum(sapply(1:m,function(i) dpois(fit_data$severe_v[fit_period,i],lambda=(model_out$sevp_v[fit_period,i]+1e-5),log=T)))
      loss.vals[6] <- sum(sapply(1:m,function(i) dpois(fit_data$severe_b[fit_period,i],lambda=(model_out$sevp_b[fit_period,i]+1e-5),log=T)))
    }
  }
  loss <- -sum(loss.vals)
  
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
  
  parms$cij_fact <- rep(1,4)
  parms$cij_fact[2:4] <- par.est[cij.fact.names]

  dr <- par.est[det.rates.names]
  parms$det_rates <- seq(dr[1],dr[2],length.out=m)
  # parms$det_rates <- c(rep(dr[1],7),rep(dr[2],2))

  if(fitUsingNB)
    parms$dp <- par.est[nbinom.dp.name]
  
  return (parms)
}

computeCij <- function(parms) {
  Cij <- lapply(1:parms$Tmax, function(t)
                mobility_data$household[t]*Cij_household*parms$cij_fact[1] +
                mobility_data$work[t]*Cij_work*parms$cij_fact[2] +
                mobility_data$school[t]*Cij_school*parms$cij_fact[3] +
                mobility_data$community[t]*Cij_community*parms$cij_fact[4])
  
  Cij <- lapply(1:length(Cij), function(t) Cij[[t]]*parms$susceptibilityFactor)
  Cij <- lapply(1:length(Cij), function(t) t(t(Cij[[t]])*parms$infectivityFactor))
  norm_fact <- parms$beta_multiplier/max(Re(eigen(Cij[[1]])$values)) #1 #
  Cij <- lapply(1:length(Cij), function(t) Cij[[t]]*norm_fact)
  return(Cij)
}

##################################################
# Fitting procedure
##################################################

if(!exists('isProfileRun'))
  isProfileRun <- FALSE

doFit <- FALSE
fitDetected <- TRUE
fitSevere <- FALSE
fitUsingNB <- TRUE


beta.mult.names <- 'R'
child.sus.factor.names <- 'delta1'
det.rates.names <- paste0('rho',1:2)
cij.fact.names <- paste0('c',1:3)
nbinom.dp.name <- 'nbinom_dp'


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

#cij factors
par.guess <- c(par.guess,rep(2,length(cij.fact.names)))
par.min <- c(par.min,rep(0,length(cij.fact.names)))
par.max <- c(par.max,rep(10,length(cij.fact.names)))
par.names <- c(par.names,cij.fact.names)

#nbinom dispersion parameter
if(fitUsingNB) {
  par.guess <- c(par.guess,rep(5,length(nbinom.dp.name)))
  par.min <- c(par.min,rep(0,length(nbinom.dp.name)))
  par.max <- c(par.max,rep(100,length(nbinom.dp.name)))
  par.names <- c(par.names,nbinom.dp.name)
}


if(!isProfileRun) {

  #best fit
  par.est <- c(2.444,0.699,0.489,0.628,1.743,1.587,2.895,5.347) #15136.7 #mean gtd=4.0, mean Pdet=3.5
  
  # best fits for sensitivity analysis
  # par.est <- c(2.332,0.718,0.527,0.714,1.259,1.283,2.239,5.457);parms$gtd <- discrete_gamma_dist(mean=3.5, sd=2.7, 14) #15111.502 #mean gtd=3.5
  # par.est <- c(2.577,0.676,0.459,0.552,2.812,2.249,4.487,5.201);parms$gtd <- discrete_gamma_dist(mean=4.5, sd=2.7, 14) #15170.422 #mean gtd=4.5
  # par.est <- c(2.425,0.701,0.502,0.643,1.674,1.491,2.877,5.335);parms$Pdet <- discrete_gamma_dist(mean=3.0, sd=2.4, 14) #15137.429 #mean Pdet=3.0
  # par.est <- c(2.470,0.694,0.474,0.611,1.847,1.723,2.925,5.350);parms$Pdet <- discrete_gamma_dist(mean=4.0, sd=2.4, 14) #15139.266 #mean Pdet=4.0
}


# par.est <- par.guess
names(par.est) <- par.names


if(isProfileRun) {
  par.est <- par.est[-profile.par.ind]
  par.min <- par.min[-profile.par.ind]
  par.max <- par.max[-profile.par.ind]
  par.names <- par.names[-profile.par.ind]
} 


if(doFit) {
  for(jj in 1:1) {
    # best <- optim(par.est, fitModel, method='Nelder-Mead',
    #               control=list(parscale=pmax(1,par.est/min(par.est[which(par.est>0)])),maxit=200,reltol=1e-4),
    #               parms=parms, fit_data=fit_data, fit_period=fit_period, par.min=par.min, par.max=par.max)
    # par.est <- best$par
    best <- optim(par.est, fitModel, method='BFGS',
    parms=parms, fit_data=fit_data, fit_period=fit_period, par.min=par.min, par.max=par.max)
    par.est <- best$par
  }
} else {
  best <- list()
  best$par <- par.est
  best$value <- fitModel(par.est, parms, fit_data, fit_period, par.min, par.max)
}

if(!isProfileRun) {
  mle_best <- best
  best.nll.val <- mle_best$value
  best.est <- mle_best$par
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
print(paste0('loss: ',round(best$value,3)))

model_out <-  run_ag_aoi_waningImmunity(parms,S0_v,S0_nv,i0_nv,i0_v,Cij,vaccinationSchedule,boostSchedule)


if(!isProfileRun) {
  
  runZ_NoBoost <-  run_ag_aoi_waningImmunity(parms,S0_v,S0_nv,i0_nv,i0_v,Cij,vaccinationSchedule,0*boostSchedule)
  
  pdead_ifsev_nv <- parms$death_rates_ifsev
  pdead_ifsev_v  <- parms$death_rates_ifsev
  pdead_ifsev_b  <- parms$death_rates_ifsev
  
  L<-parms$Tmax
  
  sim<-runZ_NoBoost;mortalityA<-sum(colSums(sim$sevp_nv[1:L,])*pdead_ifsev_nv+colSums(sim$sevp_v[1:L,])*pdead_ifsev_v+colSums(sim$sevp_b[1:L,])*pdead_ifsev_b)
  sim<-model_out;mortalityB<-sum(colSums(sim$sevp_nv[1:L,])*pdead_ifsev_nv+colSums(sim$sevp_v[1:L,])*pdead_ifsev_v+colSums(sim$sevp_b[1:L,])*pdead_ifsev_b)
  
  print(round(100*sum(runZ_NoBoost$detected[1:L,])/sum(model_out$detected[1:L,])))
  print(round(100*sum(runZ_NoBoost$sevp[1:L,])/sum(model_out$sevp[1:L,])))
  print(round(100*mortalityA/mortalityB))
}


