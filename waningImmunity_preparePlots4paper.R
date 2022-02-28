library(htmlTable)
library(magrittr)
library(kableExtra)

## This R script extract simulation data to build table 1 in the manuscript, 
## and saves releveant simulation data for use in Matlab to constructs the various plots of the paper.
## To produce the plots, use 'makeGraphs.m' in Matlab

print('Prepare plots for paper')

res <- data.frame("dates"=fitStartDate+(1:length(rowSums(fit_data$detected_nv)))-1,
                  "detected_nv"=rowSums(fit_data$detected_nv),
                  "detected_v"=rowSums(fit_data$detected_v),
                  "detected_b"=rowSums(fit_data$detected_b),
                  "severe_nv"=rowSums(fit_data$severe_nv),
                  "severe_v"=rowSums(fit_data$severe_v),
                  "severe_b"=rowSums(fit_data$severe_b))

write.csv(res,"calibrationData.csv",fileEncoding = "UTF-8")

l1 <- remove_hidden_infections_from_schedules(parms,vaccinationSchedule0,boostSchedule0)
vaccinationSchedule <- l1$vaccinationSchedule
boostSchedule <- l1$boostSchedule

l2 <- set_initial_conditions(parms,vaccinationSchedule)
i0_nv <- l2$i0_nv
i0_v  <- l2$i0_v
S0_nv <- l2$S0_nv
S0_v  <- l2$S0_v

boostScheduleEarly <- 0*boostSchedule
boostScheduleLate <- 0*boostSchedule
auxSchedule <- boostSchedule
tL <- nrow(boostSchedule)
delay <- -14 # Two weeks earlier
boostScheduleEarly[1:(tL+delay),] <- auxSchedule[(-delay+1):tL,]
delay <- 14 # Two weeks later than actual schedule
boostScheduleLate[(1+delay):tL,] <- auxSchedule[1:(tL-delay),]

### 
runZ_Boost<-  run_ag_aoi_waningImmunity(parms,S0_v,S0_nv,i0_nv,i0_v,Cij,vaccinationSchedule,boostSchedule)

print(sum(runZ_Boost$sevp))

boostScheduleEarly <- rbind(boostScheduleEarly,matrix(0,nrow=parms$Tmax,ncol=m))
boostScheduleLate <- rbind(boostScheduleLate,matrix(0,nrow=parms$Tmax,ncol=m))


### Run all simulations needed for graphs
runZ_Boost_DirectEffect <- run_ag_aoi_waningImmunity(parms,S0_v,S0_nv,i0_nv,i0_v,Cij,vaccinationSchedule,boostSchedule,isBoostProtectsAgainstTransmission=F)
runZ_Boost<-  run_ag_aoi_waningImmunity(parms,S0_v,S0_nv,i0_nv,i0_v,Cij,vaccinationSchedule,boostSchedule)
print(sum(runZ_Boost$sevp))
boostSchedule60<-boostSchedule
boostSchedule60[,1:6]<-0
runZ_Boost60 <-  run_ag_aoi_waningImmunity(parms,S0_v,S0_nv,i0_nv,i0_v,Cij,vaccinationSchedule,boostSchedule60)

boostSchedule40<-boostSchedule
boostSchedule40[,1:4]<-0
runZ_Boost40 <-  run_ag_aoi_waningImmunity(parms,S0_v,S0_nv,i0_nv,i0_v,Cij,vaccinationSchedule,boostSchedule40)

runZ_NoBoost <- run_ag_aoi_waningImmunity(parms,S0_v,S0_nv,i0_nv,i0_v,Cij,vaccinationSchedule,0*boostSchedule)

# early/late boosting
boostScheduleEarly <- 0*boostSchedule
boostScheduleLate <- 0*boostSchedule
auxSchedule <- boostSchedule
tL <- nrow(boostSchedule)
delay <- -14 # Two weeks earlier
boostScheduleEarly[1:(tL+delay),] <- auxSchedule[(-delay+1):tL,]
delay <- 14 # Two weeks later than actual schedule
boostScheduleLate[(1+delay):tL,] <- auxSchedule[1:(tL-delay),]

runZ_BoostEarly <-  run_ag_aoi_waningImmunity(parms,S0_v,S0_nv,i0_nv,i0_v,Cij,vaccinationSchedule,boostScheduleEarly)
runZ_BoostLate <-  run_ag_aoi_waningImmunity(parms,S0_v,S0_nv,i0_nv,i0_v,Cij,vaccinationSchedule,boostScheduleLate)

runZ_Boost_CI<-list()
runZ_NoBoost_CI<-list()
runZ_Boost40_CI<-list()
runZ_Boost60_CI<-list()
runZ_BoostEarly_CI<-list()
runZ_BoostLate_CI<-list()


extractDataFromSimulation <- function(sim,refSim) {

  pdead_ifsev_nv <- parms$death_rates_ifsev
  pdead_ifsev_v  <- parms$death_rates_ifsev
  pdead_ifsev_b  <- parms$death_rates_ifsev

  L<-parms$Tmax
  Lrami<-parms$Tmax-15
  detected<-sum(sim$detected[1:L,])
  severe<-sum(sim$sevp[1:L,])
  aux<-c(sum(colSums(sim$sevp_nv[1:Lrami,])*pdead_ifsev_nv),sum(colSums(sim$sevp_v[1:Lrami,])*pdead_ifsev_v),sum(colSums(sim$sevp_b[1:Lrami,])*pdead_ifsev_b))
  print(c(aux,sum(aux)))
  mortality<-sum(colSums(sim$sevp_nv[1:L,])*pdead_ifsev_nv+colSums(sim$sevp_v[1:L,])*pdead_ifsev_v+colSums(sim$sevp_b[1:L,])*pdead_ifsev_b)
  
  return(list(detected=detected,severe=severe,mortality=mortality))
}

getMinandMax <- function(res,sim) {
  minRes<-list()
  maxRes<-list()
  minRes$detected<-pmin(res$min$detected,sim$detected)
  maxRes$detected<-pmax(res$max$detected,sim$detected)
  minRes$severe<-pmin(res$min$severe,sim$severe)
  maxRes$severe<-pmax(res$max$severe,sim$severe)
  minRes$mortality<-pmin(res$min$mortality,sim$mortality)
  maxRes$mortality<-pmax(res$max$mortality,sim$mortality)
  
  return(list(min=minRes,max=maxRes))
}

res_runZ_BoostData<-list()
runZ_BoostData<-extractDataFromSimulation(runZ_Boost);
res_runZ_BoostData$min<-runZ_BoostData
res_runZ_BoostData$max<-runZ_BoostData

res_runZ_NoBoostData<-list()
runZ_NoBoostData<-extractDataFromSimulation(runZ_NoBoost,runZ_BoostData);
res_runZ_NoBoostData$min<-runZ_NoBoostData
res_runZ_NoBoostData$max<-runZ_NoBoostData

res_runZ_Boost40Data<-list()
runZ_Boost40Data<-extractDataFromSimulation(runZ_Boost40,runZ_BoostData);
res_runZ_Boost40Data$min<-runZ_Boost40Data
res_runZ_Boost40Data$max<-runZ_Boost40Data

res_runZ_Boost60Data<-list()
runZ_Boost60Data<-extractDataFromSimulation(runZ_Boost60,runZ_BoostData);
res_runZ_Boost60Data$min<-runZ_Boost60Data
res_runZ_Boost60Data$max<-runZ_Boost60Data

res_runZ_BoostEarlyData<-list()
runZ_BoostEarlyData<-extractDataFromSimulation(runZ_BoostEarly,runZ_BoostData);
res_runZ_BoostEarlyData$min<-runZ_BoostEarlyData
res_runZ_BoostEarlyData$max<-runZ_BoostEarlyData

res_runZ_BoostLateData<-list()
runZ_BoostLateData<-extractDataFromSimulation(runZ_BoostLate,runZ_BoostData);
res_runZ_BoostLateData$min<-runZ_BoostLateData
res_runZ_BoostLateData$max<-runZ_BoostLateData

if (T) {
  parameters <- readRDS('ci95_par_est.RDS')
  
for (ix in 1:length(parameters))
{
  # Prepare for run
  names(parameters[[ix]]) <- par.names
  aux_parms <- extractParameters(parameters[[ix]], parms)
  aux_Cij <- computeCij(aux_parms)

  l1 <- remove_hidden_infections_from_schedules(aux_parms,vaccinationSchedule0,boostSchedule0)
  vaccinationSchedule <- l1$vaccinationSchedule
  boostSchedule <- l1$boostSchedule
  
  aux_vaccinationSchedule <- l1$vaccinationSchedule
  aux_boostSchedule <- l1$boostSchedule
  aux_boostSchedule60<-aux_boostSchedule
  aux_boostSchedule60[,1:6]<-0
  
  aux_boostSchedule40<-aux_boostSchedule
  aux_boostSchedule40[,1:4]<-0
  
  aux_boostScheduleEarly <- 0*boostSchedule
  aux_boostScheduleLate <- 0*boostSchedule
  tL <- nrow(boostSchedule)
  
  delay <- -14 # Two weeks earlier
  boostScheduleEarly[1:(tL+delay),] <- auxSchedule[(-delay+1):tL,]
  delay <- 14 # Two weeks later than actual schedule
  boostScheduleLate[(1+delay):tL,] <- auxSchedule[1:(tL-delay),]
  
  l2 <- set_initial_conditions(aux_parms,aux_vaccinationSchedule)
  i0_nv <- l2$i0_nv
  i0_v  <- l2$i0_v
  S0_nv <- l2$S0_nv
  S0_v  <- l2$S0_v
  
  ### Run simulations
  res_runZ_Boost_CI <-  extractDataFromSimulation(run_ag_aoi_waningImmunity(aux_parms,S0_v,S0_nv,i0_nv,i0_v,aux_Cij,aux_vaccinationSchedule,aux_boostSchedule))
  res_runZ_NoBoost_CI <-  extractDataFromSimulation(run_ag_aoi_waningImmunity(aux_parms,S0_v,S0_nv,i0_nv,i0_v,aux_Cij,aux_vaccinationSchedule,0*boostSchedule),res_runZ_Boost_CI)
  res_runZ_Boost40_CI <-  extractDataFromSimulation(run_ag_aoi_waningImmunity(aux_parms,S0_v,S0_nv,i0_nv,i0_v,aux_Cij,aux_vaccinationSchedule,aux_boostSchedule40),res_runZ_Boost_CI)
  res_runZ_Boost60_CI <-  extractDataFromSimulation(run_ag_aoi_waningImmunity(aux_parms,S0_v,S0_nv,i0_nv,i0_v,aux_Cij,aux_vaccinationSchedule,aux_boostSchedule60),res_runZ_Boost_CI)
  res_runZ_BoostEarly_CI <-  extractDataFromSimulation(run_ag_aoi_waningImmunity(aux_parms,S0_v,S0_nv,i0_nv,i0_v,aux_Cij,aux_vaccinationSchedule,boostScheduleEarly),res_runZ_Boost_CI)
  res_runZ_BoostLate_CI<-extractDataFromSimulation(run_ag_aoi_waningImmunity(aux_parms,S0_v,S0_nv,i0_nv,i0_v,aux_Cij,aux_vaccinationSchedule,boostScheduleLate),res_runZ_Boost_CI)

  res_runZ_BoostData<-getMinandMax(res_runZ_BoostData,res_runZ_Boost_CI)
  res_runZ_NoBoostData<-getMinandMax(res_runZ_NoBoostData,res_runZ_NoBoost_CI)
  res_runZ_Boost40Data<-getMinandMax(res_runZ_Boost40Data,res_runZ_Boost40_CI)
  res_runZ_Boost60Data<-getMinandMax(res_runZ_Boost60Data,res_runZ_Boost60_CI)
  res_runZ_BoostEarlyData<-getMinandMax(res_runZ_BoostEarlyData,res_runZ_BoostEarly_CI)
  res_runZ_BoostLateData<-getMinandMax(res_runZ_BoostLateData,res_runZ_BoostLate_CI)
  
  print(c(res_runZ_BoostData$min$severe,res_runZ_Boost_CI$severe,res_runZ_BoostData$max$severe))
}

save.image(file ='SimulationData4Paper.Rdata')

load('SimulationData4Paper.Rdata')

### Table presenting reduction of epidemiological outcomes in various scenarios
tab<-rep(0,6)
tabMin<-rep(0,6)
tabMax<-rep(0,6)
ix<-1
L<-150
a<-1;ix<-1
tab[ix]<-1-runZ_NoBoostData$detected/a;tabMin[ix]<-1-res_runZ_NoBoostData$min$detected/a;tabMax[ix]<-1-res_runZ_NoBoostData$max$detected/a;ix<-ix+1
tab[ix]<-1-runZ_Boost60Data$detected/a;tabMin[ix]<-1-res_runZ_Boost60Data$min$detected/a;tabMax[ix]<-1-res_runZ_Boost60Data$max$detected/a;ix<-ix+1
tab[ix]<-1-runZ_Boost40Data$detected/a;tabMin[ix]<-1-res_runZ_Boost40Data$min$detected/a;tabMax[ix]<-1-res_runZ_Boost40Data$max$detected/a;ix<-ix+1
tab[ix]<-1-runZ_BoostData$detected/a;tabMin[ix]<-1-res_runZ_BoostData$min$detected/a;tabMax[ix]<-1-res_runZ_BoostData$max$detected/a;ix<-ix+1
tab[ix]<-1-runZ_BoostEarlyData$detected/a;tabMin[ix]<-1-res_runZ_BoostEarlyData$min$detected/a;tabMax[ix]<-1-res_runZ_BoostEarlyData$max$detected/a;ix<-ix+1
tab[ix]<-1-runZ_BoostLateData$detected/a;tabMin[ix]<-1-res_runZ_BoostLateData$min$detected/a;tabMax[ix]<-1-res_runZ_BoostLateData$max$detected/a;ix<-ix+1
detected<-1-tab
detected_min<-1-tabMin
detected_max<-1-tabMax

a<-1;ix<-1
tab[ix]<-1-runZ_NoBoostData$severe/a;tabMin[ix]<-1-res_runZ_NoBoostData$min$severe/a;tabMax[ix]<-1-res_runZ_NoBoostData$max$severe/a;ix<-ix+1
tab[ix]<-1-runZ_Boost60Data$severe/a;tabMin[ix]<-1-res_runZ_Boost60Data$min$severe/a;tabMax[ix]<-1-res_runZ_Boost60Data$max$severe/a;ix<-ix+1
tab[ix]<-1-runZ_Boost40Data$severe/a;tabMin[ix]<-1-res_runZ_Boost40Data$min$severe/a;tabMax[ix]<-1-res_runZ_Boost40Data$max$severe/a;ix<-ix+1
tab[ix]<-1-runZ_BoostData$severe/a;tabMin[ix]<-1-res_runZ_BoostData$min$severe/a;tabMax[ix]<-1-res_runZ_BoostData$max$severe/a;ix<-ix+1
tab[ix]<-1-runZ_BoostEarlyData$severe/a;tabMin[ix]<-1-res_runZ_BoostEarlyData$min$severe/a;tabMax[ix]<-1-res_runZ_BoostEarlyData$max$severe/a;ix<-ix+1
tab[ix]<-1-runZ_BoostLateData$severe/a;tabMin[ix]<-1-res_runZ_BoostLateData$min$severe/a;tabMax[ix]<-1-res_runZ_BoostLateData$max$severe/a;ix<-ix+1
severe <- 1-tab
severe_min <- 1-tabMin
severe_max <- 1-tabMax

a<-1;ix<-1
tab[ix]<-1-runZ_NoBoostData$mortality/a;tabMin[ix]<-1-res_runZ_NoBoostData$min$mortality/a;tabMax[ix]<-1-res_runZ_NoBoostData$max$mortality/a;ix<-ix+1
tab[ix]<-1-runZ_Boost60Data$mortality/a;tabMin[ix]<-1-res_runZ_Boost60Data$min$mortality/a;tabMax[ix]<-1-res_runZ_Boost60Data$max$mortality/a;ix<-ix+1
tab[ix]<-1-runZ_Boost40Data$mortality/a;tabMin[ix]<-1-res_runZ_Boost40Data$min$mortality/a;tabMax[ix]<-1-res_runZ_Boost40Data$max$mortality/a;ix<-ix+1
tab[ix]<-1-runZ_BoostData$mortality/a;tabMin[ix]<-1-res_runZ_BoostData$min$mortality/a;tabMax[ix]<-1-res_runZ_BoostData$max$mortality/a;ix<-ix+1
tab[ix]<-1-runZ_BoostEarlyData$mortality/a;tabMin[ix]<-1-res_runZ_BoostEarlyData$min$mortality/a;tabMax[ix]<-1-res_runZ_BoostEarlyData$max$mortality/a;ix<-ix+1
tab[ix]<-1-runZ_BoostLateData$mortality/a;tabMin[ix]<-1-res_runZ_BoostLateData$min$mortality/a;tabMax[ix]<-1-res_runZ_BoostLateData$max$mortality/a;ix<-ix+1
mortality <- 1-tab
mortality_min <- 1-tabMin
mortality_max <- 1-tabMax

rowNames <- c("No Boost","Ages 60 and above","Ages 40 and above","Actual","Early schedule","Late schedule")
colNames <- c("Booster Policy","Change in cases","Change in severe cases","Change in mortality")

db <- data.frame(rowNames,round(detected_min),round(detected), round(detected_max),round(1*severe_min),round(1*severe), round(1*severe_max),round(1*mortality_min),round(1*mortality), round(1*mortality_max))
kbl(db,"latex")
kbl(db)
  
refDetected=runZ_BoostData$detected
refSevere=runZ_BoostData$severe
refMortality=runZ_BoostData$mortality
db <- data.frame(rowNames,paste(round(100*detected/refDetected),  '% (',round(100*detected_min/refDetected),  '-',round(100*detected_max/refDetected),')',sep=''),
                          paste(round(100*severe/refSevere),      '% (',round(100*severe_min/refSevere),      '-',round(100*severe_max/refSevere),')',sep=''),
                          paste(round(100*mortality/refMortality),'% (',round(100*mortality_min/refMortality),'-',round(100*mortality_max/refMortality),')',sep=''))
#db <- data.frame(rowNames,paste((100*detected/refDetected),  '% (',(100*detected_min/refDetected),  '-',(100*detected_max/refDetected),')',sep=''),
#                 paste((100*severe/refSevere),      '% (',(100*severe_min/refSevere),      '-',(100*severe_max/refSevere),')',sep=''),
#                 paste((100*mortality/refMortality),'% (',(100*mortality_min/refMortality),'-',(100*mortality_max/refMortality),')',sep=''))

kbl(db,col.names=colNames,"latex")

}
### CalibrationGraph 

t=simulationStartDate-1+1:parms$Tmax

res <- data.frame("dates"=t,
                  "detected_nv"=rowSums(runZ_Boost$detected_nv),
                  "detected_v"=rowSums(runZ_Boost$detected_v),
                  "detected_b"=rowSums(runZ_Boost$detected_b),
                  "severe_nv"=rowSums(runZ_Boost$sevp_nv),
                  "severe_v"=rowSums(runZ_Boost$sevp_v),
                  "severe_b"=rowSums(runZ_Boost$sevp_b))
write.csv(res,"calibration.csv",fileEncoding = "UTF-8")

res <- data.frame("dates"=t,
                  "detected_NoBoost"=rowSums(runZ_NoBoost$detected),
                  "detected_Boost"=rowSums(runZ_Boost$detected),
                  "detected_Boost_DirectEffect_overall"=rowSums(runZ_Boost_DirectEffect$detected),
                  "detected_Boost_DirectEffect_v3"=rowSums(runZ_Boost_DirectEffect$detected_b),
                  "detected_Boost_v3"=rowSums(runZ_Boost$detected_b),
                  "detected_Boost60"=rowSums(runZ_Boost60$detected),
                  "detected_Boost40"=rowSums(runZ_Boost40$detected),
                  "severe_NoBoost"=rowSums(runZ_NoBoost$sevp),
                  "severe_Boost_DirectEffect_overall"=rowSums(runZ_Boost_DirectEffect$sevp),
                  "severe_Boost_DirectEffect_v3"=rowSums(runZ_Boost_DirectEffect$sevp_b),
                  "severe_Boost"=rowSums(runZ_Boost$sevp),
                  "severe_Boost_v3"=rowSums(runZ_Boost$sevp_b),
                  "severe_Boost60"=rowSums(runZ_Boost60$sevp),
                  "severe_Boost40"=rowSums(runZ_Boost40$sevp),
                  "R_NoBoost"=runZ_NoBoost$Reff_observed,
                  "R_Boost"=runZ_Boost$Reff_observed,
                  "R_Boost60"=runZ_Boost60$Reff_observed,
                  "R_Boost40"=runZ_Boost40$Reff_observed,
                  "BoostData60"=rowSums(boostSchedule[(dv+1):(dv+parms$Tmax),7:9]),
                  "BoostData4059"=rowSums(boostSchedule[(dv+1):(dv+parms$Tmax),5:6]),
                  "BoostData039"=rowSums(boostSchedule[(dv+1):(dv+parms$Tmax),1:4])
                  )
write.csv(res,"EffectOfBoost.csv",fileEncoding = "UTF-8")


res <- data.frame("dates"=t,
                  "Reff_observed"=runZ_Boost$Reff_observed,
                  "Reff"=runZ_Boost$Reff,
                  "R_early"=runZ_BoostEarly$Reff_observed,
                  "R_late_observed"=runZ_BoostLate$Reff_observed,
                  "R_late"=runZ_BoostLate$Reff_observed,
                  "detected"=rowSums(runZ_Boost$detected),
                  "detected_Early"=rowSums(runZ_BoostEarly$detected),
                  "detected_Late"=rowSums(runZ_BoostLate$detected),
                  "severe"=rowSums(runZ_Boost$sevp),
                  "severe_Early"=rowSums(runZ_BoostEarly$sevp),
                  "severe_Late"=rowSums(runZ_BoostLate$sevp)
                  )
write.csv(res,"longRangeProjection.csv",fileEncoding = "UTF-8")


res <- data.frame("dates"=t,
                  "detected"=rowSums(runZ_Boost$detected),
                  "detected_Early"=rowSums(runZ_BoostEarly$detected),
                  "detected_Late"=rowSums(runZ_BoostLate$detected),
                  "severe"=rowSums(runZ_Boost$sevp),
                  "severe_Early"=rowSums(runZ_BoostEarly$sevp),
                  "severe_Late"=rowSums(runZ_BoostLate$sevp),
                  "R_early"=runZ_BoostEarly$Reff_observed,
                  "Reff"=runZ_Boost$Reff_observed,
                  "R_late"=runZ_BoostLate$Reff_observed
)
write.csv(res,"EffectOfTiming.csv",fileEncoding = "UTF-8")

print('Finished preparing plots for paper')
