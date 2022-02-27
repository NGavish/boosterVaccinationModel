

discrete_gamma_dist <- function(mean, sd, d) {
  x <- 1:d
  shape <- (mean/sd)^2
  scale <- mean/shape
  P <- pgamma(x+1,shape=shape,scale=scale)-pgamma(x,shape=shape,scale=scale)
  P <- P/sum(P)
  return(P)
}

gather_parameters <- function() {
  
  ######################
  # Gather parameters
  ######################

  ## Demography
  ag_labels <- c('0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+')
  N <- c(1830000,1550000,1310000,1210000,1106000,864000,747000,500000,283000)

  ag_dist <- N/sum(N)
  m <- length(N)
  
  gtd <- discrete_gamma_dist(mean=4.0, sd=2.7, 14)   # generation time distribution
  Pdet <- discrete_gamma_dist(mean=3.5, sd=2.4, 14)  # detection time distribution
  
  # initial detection rates - to be fitted
  det_rates <- rep(0.5,m)

  death_rates_ifsev <- c(0.02, 0.00, 0.05, 0.07, 0.08, 0.14, 0.30, 0.37, 0.50) 
  
  # (age-dependent susceptibility and infectivity)
  infectivityFactor <- rep(1,m)   
  susceptibilityFactor <- rep(1,m)  
  susceptibilityFactor[1] <- 0.5
  
  Tv1 <- 7             #Duration of phase 0
  Tv2 <- 28-Tv1        #Duration of phase 1
  Tv3 <- 365-Tv2-Tv1   #Duration of phase 3
  dv <- Tv1+Tv2+Tv3
  
  #VE Profile
  SeffTemporalprofile <- rep(0,dv)
  SeffTemporalprofile[1:34] <- seq(0,0.97,length.out=34)
  SeffTemporalprofile[34+(1:300)] <- seq(0.97,0,length.out=300)
  effInfProfile <- SeffTemporalprofile%*%t(rep(1,m))
  
  # Boost decays as vaccine
  bSeffTemporalprofile <- rep(0,dv)
  VE0 <- SeffTemporalprofile[200]
  factor <- c(rep(2,6),seq(2,15),rep(15,14))
  bSeffTemporalprofile[1:34]<-1-(1-VE0)/factor
  bSeffTemporalprofile[35:dv] <- bSeffTemporalprofile[34]*SeffTemporalprofile[35:dv]/SeffTemporalprofile[35]

  # slow wanning - Boost decays half as fast as vaccine
  bSeffTemporalprofileSlowWaning <- bSeffTemporalprofile
  TL <- length(bSeffTemporalprofile)-35
  bSeffTemporalprofileSlowWaning[35+0:TL] <- approx(35+0:TL,bSeffTemporalprofile[35+0:TL],35+(0:TL)/2)$y
  beffInfProfileSlowWaning <- bSeffTemporalprofileSlowWaning%*%t(rep(1,m))
  
  # known unvaccinated and vaccinated infected on July 1st 2021 by age group and time since vaccination
  i0_known_nv <- readRDS('i0_nv.RDS')
  i0_known_v <- readRDS('i0_v.RDS')
  
  # unvaccinated and vaccinated recovered up to July 1st 2021 by age group
  recovered_known_nv <- c(110906,172612,161362,116646,97847,73210,47565,24320,16793)
  recovered_known_v <-  c(3,59,356,654,931,815,909,746,648)

  return(list(N=N,m=m,ag_dist=ag_dist,ag_labels=ag_labels,gtd=gtd,Pdet=Pdet,
          det_rates=det_rates,death_rates_ifsev=death_rates_ifsev,
          susceptibilityFactor=susceptibilityFactor,infectivityFactor=infectivityFactor,
          Tv1=Tv1,Tv2=Tv2,Tv3=Tv3,dv=dv,
          effInfProfile=effInfProfile, beffInfProfile=beffInfProfileSlowWaning,
  			  i0_known_nv=i0_known_nv,i0_known_v=i0_known_v,
  			  recovered_known_nv=recovered_known_nv,recovered_known_v=recovered_known_v))
}

smooth_window <- function(dat, win_size) {
  zoo::rollapply(dat,FUN=mean,width=win_size,align="center")
}

#######################################################
# Parameters, initial conditions and other definitions
#######################################################

# Gather parameters & extract out common parameters 
parms <- gather_parameters();
m     <- parms$m # Number of age groups
N     <- parms$N # Size of age groups
dv    <- parms$dv
d     <- length(parms$gtd) 

simulationStartDate <- as.Date("2021-07-01")
simulationEndDate <- as.Date("2021-12-31")

Tmax <- as.numeric(simulationEndDate-simulationStartDate)
parms$Tmax <- Tmax

dates <-  simulationStartDate+(0:Tmax)

##################################################
# Read detailed data used for fitting
##################################################


detectedData <- readRDS('detected.RDS')
severeData <- readRDS('severe.RDS')

fit_data <- list()
fit_data$detected_nv <- detectedData$nv
fit_data$detected_v <- detectedData$v
fit_data$detected_b <- detectedData$b
fit_data$severe_nv <- severeData$nv
fit_data$severe_v <- severeData$v
fit_data$severe_b <- severeData$b

fitStartDate <- simulationStartDate 
fitEndDate <- as.Date("2021-11-25") 

fit_length <- as.numeric(fitEndDate-fitStartDate+1)
fit_period <- 1:fit_length
fit_dates <- fitStartDate+1:fit_length-1


parms$sevp_rates_nv <- colSums(fit_data$severe_nv)/colSums(fit_data$detected_nv)
parms$sevp_rates_v <- colSums(fit_data$severe_v)/colSums(fit_data$detected_v)
parms$sevp_rates_b <- colSums(fit_data$severe_b)/colSums(fit_data$detected_b)
parms$sevp_rates_b[1] <- 0

##############################################
# Mobility data
##############################################

mobility_data0 <- read.csv('Mobility_Report_Google.csv')
mobility_data <- data.frame(household=(1+mobility_data0$residential_percent_change_from_baseline/100),
                            work=(1+mobility_data0$workplaces_percent_change_from_baseline/100),
                            community=(1+mobility_data0$retail_and_recreation_percent_change_from_baseline/100),
                            school=(1+mobility_data0$schools/100))

mobility_data$household[4:(length(mobility_data$household)-3)] <- smooth_window(mobility_data$household,7)
mobility_data$work[4:(length(mobility_data$work)-3)] <- smooth_window(mobility_data$work,7)
mobility_data$school[4:(length(mobility_data$school)-3)] <- smooth_window(mobility_data$school,7)
mobility_data$community[4:(length(mobility_data$community)-3)] <- smooth_window(mobility_data$community,7)

if(parms$Tmax>nrow(mobility_data)) {
  mobility_data[(nrow(mobility_data)+1):(parms$Tmax),] <- mobility_data[nrow(mobility_data),]
}

##############################################
# Vaccination & boost schedule
##############################################

# vaccinationSchedule contains detailed vaccination data up to dv days prior to the simulation start time

vaccinationTable <- readRDS('vaccinationTable.RDS')
vacDataStartDate <- as.Date(vaccinationTable[1,"VaccinationDate"][[1]],format='%Y-%m-%d')
vacDataEndDate <- as.Date(vaccinationTable[nrow(vaccinationTable),"VaccinationDate"][[1]],format='%Y-%m-%d')
vaccinationSchedule <- as.matrix(vaccinationTable[,2:10])

boostTable <- readRDS('boostTable.RDS')
boostSchedule <- as.matrix(boostTable[,2:10])

# If relevant, aggregate all detailed vaccination data more than dv days prior to the simulation start date to last day detailed records are kept
if ((simulationStartDate-vacDataStartDate)>dv) { 
  idx=as.numeric(simulationStartDate-vacDataStartDate-dv)
  vaccinationSchedule[idx,]<-colSums(vaccinationSchedule[1:idx,])
  vaccinationSchedule<-vaccinationSchedule[idx:nrow(vaccinationSchedule),]
} else # Otherwise, add trailing zeros so that vaccinationSchedule accounts for dv days prior to simulation start date 
{
  th<-as.numeric(vacDataStartDate-simulationStartDate+dv-1)
  vaccinationSchedule<-rbind(matrix(0,th,m),vaccinationSchedule)
  boostSchedule<-rbind(matrix(0,th,m),boostSchedule)
}

# If relevant, pad the vaccination schedule with zeros until the end of the simulation
if (vacDataEndDate<(simulationStartDate+Tmax)) {
  vaccinationSchedule<-rbind(vaccinationSchedule,matrix(0,simulationStartDate+Tmax-vacDataEndDate+1,m))
  boostSchedule<-rbind(boostSchedule,matrix(0,simulationStartDate+Tmax-vacDataEndDate+1,m))
}

remove_hidden_infections_from_schedules <- function(parms,vaccinationSchedule,boostSchedule) {
  # Estimate numbers of unknown recovered per age group, 
  # and reduce unknown recovered that were vaccinated from the vaccination schedule
  recovered_known   <- parms$recovered_known_nv+parms$recovered_known_v
  recovered_overall <- recovered_known/parms$det_rates
  recovered_unknown <- recovered_overall-recovered_known
  hiddenPercent     <- recovered_unknown/(N-recovered_known)
  vaccinationSchedule <- t(t(vaccinationSchedule)*(1-hiddenPercent))
  boostSchedule <- t(t(boostSchedule)*(1-hiddenPercent))
  return (list(vaccinationSchedule=vaccinationSchedule,boostSchedule=boostSchedule))
}

vaccinationSchedule0 <- vaccinationSchedule
boostSchedule0 <- boostSchedule
l1 <- remove_hidden_infections_from_schedules(parms,vaccinationSchedule0,boostSchedule0)
vaccinationSchedule <- l1$vaccinationSchedule
boostSchedule <- l1$boostSchedule


##################################################
# Initial Conditions
##################################################

set_initial_conditions <- function(parms,vaccinationSchedule) {
  
  d <- length(parms$gtd)
  recovered_known   <- parms$recovered_known_nv+parms$recovered_known_v
  recovered_overall <- recovered_known/parms$det_rates 
  recovered_unknown <- recovered_overall-recovered_known
  
  i0_nv <- parms$i0_known_nv/parms$det_rates
  i0_v <- parms$i0_known_v/parms$det_rates
  
  S0_v <- matrix(0,dv,m)
  S0_v[1:(dv-1),] <- vaccinationSchedule[(dv-1):1,]
  
  # Reduce those that were infected since their vaccination (already removed the unknown infected from vaccinationSchedule)
  vac_not_inf <- pmax(0,1-parms$recovered_known_v/colSums(vaccinationSchedule[1:(dv-1),]))
  S0_v <- sapply(1:m,function(i) S0_v[,i]*vac_not_inf[i]) 
  S0_nv <- parms$N-recovered_overall-colSums(i0_nv)-colSums(vaccinationSchedule[1:(dv-1),])-colSums(apply(i0_v,2,rowSums))
  
  return (list(i0_nv=i0_nv,i0_v=i0_v,S0_nv=S0_nv,S0_v=S0_v))
}

l2 <- set_initial_conditions(parms,vaccinationSchedule)
i0_nv <- l2$i0_nv
i0_v  <- l2$i0_v
S0_nv <- l2$S0_nv
S0_v  <- l2$S0_v

##################################################
# Contact matrices
##################################################

Cij_household <- readRDS('./contact_mats/mat_household.RDS')
Cij_work      <- readRDS('./contact_mats/mat_work.RDS')
Cij_school    <- readRDS('./contact_mats/mat_school.RDS')
Cij_community <- readRDS('./contact_mats/mat_community.RDS')


