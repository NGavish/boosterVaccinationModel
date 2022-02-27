
run_ag_aoi_waningImmunity <- function(parms,  # General set of parameters
                                      S0_v,S0_nv,i0_nv,i0_v, # Initial conditions
                                      Cij, # Contact matrix
                                      vaccinationSchedule,
                                      boostSchedule,
                                      isBoostProtectsAgainstTransmission=T) {
  
  Tmax<- parms$Tmax   # Simulation end time
  N   <- parms$N      # Size of age groups
  m   <- parms$m      # Number of age groups
  gtd <- parms$gtd    # generation time distribution
  
  #detection rates
  det_rates_nv <- parms$det_rates
  det_rates_v <- parms$det_rates
  det_rates_b <- parms$det_rates
  
  #severe rates
  sevp_rates_nv<-parms$sevp_rates_nv
  sevp_rates_v<-parms$sevp_rates_v
  sevp_rates_b<-parms$sevp_rates_b
  
  Tv1 <- parms$Tv1
  Tv2 <- parms$Tv2
  Tv3 <- parms$Tv3
  dv <- parms$dv
  
  d    <- length(parms$gtd) 
  ddet <- length(parms$Pdet)

  # Efficacy profiles as function of time since vaccine/boost and age group
  effInfProfile <- parms$effInfProfile
  beffInfProfile <- parms$beffInfProfile
  if(!isBoostProtectsAgainstTransmission)
    beffInfProfile <- effInfProfile
  
  # Compartments
  infected_nv <- matrix(0,Tmax+dv,m)
  infected_v <- array(0,c(Tmax+dv,m,dv)) # (t,j,s) - Vaccinated at day s and infected on day t
  infected_b <- array(0,c(Tmax+dv,m,dv)) # (t,j,s) - Got a boost at day s and infected on day t
  
  sevp <- matrix(0,Tmax+dv,m)
  sevp_nv <- matrix(0,Tmax+dv,m)
  sevp_v <- matrix(0,Tmax+dv,m)   
  sevp_b <- matrix(0,Tmax+dv,m)
  
  detected <- matrix(0,Tmax+dv,m)
  detected_nv <- matrix(0,Tmax+dv,m)
  detected_v <- array(0,c(Tmax+dv,m))
  detected_b <- array(0,c(Tmax+dv,m))
  
  sus_nv <- rep(0,m)      # [j]   - Overall number of susceptible that had not been vaccinated 
  sus_v <- matrix(0,m,dv) # [j,s] - Vaccinated at day s and had not been infected 
  sus_b <- matrix(0,m,dv) # [j,s] - Got a boost on day s and had not been infected 
  
  # Initializations
  infected_nv[(1+dv-d):dv,] <- i0_nv  # [t,j] - Daily number of newly infected at day t that had not been vaccinated 
  infected_v[(1+dv-d):dv,,] <- i0_v   # [t,j,s] - Daily number of newly infected at day t, who had been vaccinated on day s
  
  sus_nv <- S0_nv
  sus_v[,1:nrow(S0_v)] <- t(S0_v)
  
  # Additional variables
  Reff <- rep(0,Tmax)
  Reff_observed <- rep(0,Tmax)
  
  for(t in (dv+1):(Tmax+dv)) {
    
    infected <- infected_nv[(t-d):(t-1),] +
      apply(infected_v[(t-d):(t-1),,],2,rowSums) +
      apply(infected_b[(t-d):(t-1),,],2,rowSums)
    
    finf <- rowSums(sapply(1:m, function(k) Cij[[t-dv]][,k]/N[k]*sum(gtd[d:1]*infected[,k])))
    
    infected_nv[t,] <- pmax(0,sus_nv*finf)
    infected_v[t,,] <- pmax(0,t((1-effInfProfile)*t(sus_v*finf)))
    infected_b[t,,] <- pmax(0,t((1-beffInfProfile)*t(sus_b*finf)))
    
    # Detected & Severe 
    t1.det <- max(1,(t-ddet+1)):t
    detected_nv[t,] <- rowSums(det_rates_nv*t(parms$Pdet[length(t1.det):1]*infected_nv[t1.det,]))
    for (j in 1:m) {
      detected_nv[t,j] <- detected_nv[t,j]+sum(det_rates_nv[j]*parms$Pdet[length(t1.det):1]*rowSums(infected_v[t1.det,j,1:(Tv1+Tv2-1)])) # Recently vaccinated are counted as non-vaccinated
      detected_v[t,j] <- sum(det_rates_v[j]*parms$Pdet[length(t1.det):1]*rowSums(infected_v[t1.det,j,(Tv1+Tv2):dv]))
      detected_b[t,j] <- sum(det_rates_b[j]*parms$Pdet[length(t1.det):1]*rowSums(infected_b[t1.det,j,1:dv])) 
    }
    sevp_nv[t,] <- (sevp_rates_nv*detected_nv[t,])
    sevp_v[t,] <- (sevp_rates_v*detected_v[t,])
    sevp_b[t,] <- (sevp_rates_b*detected_b[t,])
    
    #Reff calculations
    detected[t,] <- detected_nv[t,]+detected_v[t,]+detected_b[t,]
    sus_aux <- (sus_nv+rowSums(t((1-effInfProfile))*sus_v)+rowSums(t((1-beffInfProfile))*sus_b))/N
    Reff[t-dv] <- max(Re(eigen(t(Cij[[t-dv]])*sus_aux)$values))
    Reff_observed[t-dv] <- sum(detected[t,])/sum(sapply(1:m, function(k) sum(gtd[d:1]*(detected[(t-d):(t-1),k]))))
    
    # Vaccination & book-keeping for the next day
    vac_t <- pmin(vaccinationSchedule[t,],sus_nv-infected_nv[t,])
    sus_nv <- sus_nv-infected_nv[t,]-vac_t
    sus_v <- sus_v - infected_v[t,,]
    sus_b <- sus_b - infected_b[t,,]
    
    aux <- sus_v[,dv] 
    sus_v[,2:dv] <- sus_v[,1:(dv-1)]
    sus_v[,dv]<-sus_v[,dv]+aux 
    sus_v[,1] <- vac_t 
    
    aux <- sus_b[,dv] 
    sus_b[,2:dv] <- sus_b[,1:(dv-1)]
    sus_b[,dv]<-sus_b[,dv]+aux 
    sus_b[,1]<-0 # Code for administering boost is next
    
    # Administer boost
    for (j in 1:m) {
      # First give boost to those fully vaccinated
      boostLeft <- boostSchedule[t,j]
      # Assume FIFO - most optimistic case
      st<-dv
      while ((boostLeft>0) & (st>0))
      {
        boostGiven <- min(sus_v[j,st],boostLeft)
        sus_v[j,st]<-sus_v[j,st]-boostGiven
        if(isBoostProtectsAgainstTransmission)
          sus_b[j,1]<-sus_b[j,1]+boostGiven
        else
          sus_b[j,st]<-sus_b[j,st]+boostGiven
        boostLeft <- boostLeft-boostGiven
        st <- st-1
      }
    }
  }
  sus_v <- rowSums(sus_v)
  sus_b <- rowSums(sus_b)
  
  infected_v <- apply(infected_v,2,rowSums)
  infected_b <- apply(infected_b,2,rowSums)
  
  inf_nv <- as.matrix(infected_nv[(dv+1):(Tmax+dv),])
  inf_v <- as.matrix(infected_v[(dv+1):(Tmax+dv),])
  inf_b <- as.matrix(infected_b[(dv+1):(Tmax+dv),])
  
  inf <- inf_nv+inf_v+inf_b
  
  det_nv <- as.matrix(detected_nv[(dv+1):(Tmax+dv),])
  det_v <- as.matrix(detected_v[(dv+1):(Tmax+dv),])
  det_b <- as.matrix(detected_b[(dv+1):(Tmax+dv),])
  det <- det_nv+det_v+det_b
  
  sevp_nv <- as.matrix(sevp_nv[(dv+1):(Tmax+dv),])
  sevp_v <- as.matrix(sevp_v[(dv+1):(Tmax+dv),])
  sevp_b <- as.matrix(sevp_b[(dv+1):(Tmax+dv),])
  sevp <- sevp_nv+sevp_v+sevp_b
  
  return (list(sus_nv=sus_nv,sus_v=sus_v,sus_b=sus_b,
               infected_nv=inf_nv,infected_v=inf_v,infected=inf,
               detected_nv=det_nv,detected_v=det_v,detected_b=det_b,detected=det,
               sevp_nv=sevp_nv,sevp_v=sevp_v,sevp_b=sevp_b,sevp=sevp,
               Reff=Reff,Reff_observed=Reff_observed))
}