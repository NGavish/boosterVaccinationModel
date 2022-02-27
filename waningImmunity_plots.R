library(fields)
library(scales)

x11(width=8);par(mfcol=c(3,2))
barplot(colSums(fit_data$detected_nv),names.arg=parms$ag_labels,ylab='detected',main='unvaccinated',cex.axis=1.2,cex.lab=1.2)
barplot(colSums(fit_data$detected_v),names.arg=parms$ag_labels,ylab='detected',main='vaccinated',cex.axis=1.2,cex.lab=1.2)
barplot(colSums(fit_data$detected_b),names.arg=parms$ag_labels,ylab='detected',main='boostered',cex.axis=1.2,cex.lab=1.2)
barplot(colSums(fit_data$severe_nv),names.arg=parms$ag_labels,ylab='severe',main='unvaccinated',cex.axis=1.2,cex.lab=1.2)
barplot(colSums(fit_data$severe_v),names.arg=parms$ag_labels,ylab='severe',main='vaccinated',cex.axis=1.2,cex.lab=1.2)
barplot(colSums(fit_data$severe_b),names.arg=parms$ag_labels,ylab='severe',main='boostered',cex.axis=1.2,cex.lab=1.2)

# x11();
# barplot(rbind(parms$sevp_rates_nv,parms$sevp_rates_v,parms$sevp_rates_b),
#         beside=T,col=2:4,names.arg=parms$ag_labels,ylim=c(0,0.5),ylab='probability of detected cases developing severe illness')
# grid(nx=NA,ny=NULL)
# legend('topleft',legend=c('unvaccinated','vaccinated','boostered'),fill=2:4)


x11();image.plot(Cij[[1]], axes=F, useRaster=T, xlab='Age of Individual',ylab='Age of Contact')

x11(width=10); matplot(mobility_data[fit_period,],type='l',lty=1:4,lwd=3,xaxt='n',cex.axis=1.5,ylab='1-(% change from baseline)/100 ',cex.axis=1.5,cex.lab=1.5)
axis(1, at=seq(1,fit_length,by=14), labels=format(dates[seq(1,fit_length,by=14)],'%d/%m'), las=2 ,cex.axis=1.5)
legend('bottomright',c('household','work','school','community'),lwd=rep(3,4),lty=c(1,2,4,3),col=c(1,2,4,3),cex=1.5)

x11(); plot(model_out$Reff,col='blue',lwd=3,type='l',xaxt='n',xlab='',ylab='R',ylim=c(0.7,1.5),cex.axis=1.5,cex.lab=1.5);
lines(14:length(model_out$Reff_observed),model_out$Reff_observed[14:length(model_out$Reff_observed)],lwd=2)
axis(1, at=seq(1,length(model_out$Reff_observed),by=7), labels=format(dates[seq(1,length(model_out$Reff_observed),by=7)],'%d/%m'), las=2 ,cex.axis=1.5)
grid(nx=NA,ny=NULL)


# Daily detected and severe per vac status and age-group
x11(height=10); par(mfrow=c(3,3))
for(i in 1:m) {
  plot(fit_data$detected_nv[,i],xaxt='n',xlab='',ylab='new detected cases',main=paste0('unvaccinated (',parms$ag_labels[i],')'),cex.main=1.25,cex.axis=1.25,cex.lab=1.25,
       ylim=c(0,max(model_out$detected_nv_high))); 
  polygon(c(1:Tmax,Tmax:1),c(model_out$detected_nv_high[,i],rev(model_out$detected_nv_low[,i])),col='lightgray')
  points(fit_data$detected_nv[,i])
  lines(model_out$detected_nv[fit_period,i],col='red',lwd=3)
  axis(1, at=seq(1,nrow(fit_data$detected_nv),by=7), labels=format(dates[seq(1,nrow(fit_data$detected_nv),by=7)],'%d/%m'), las=2 ,cex.axis=1.25)
}
x11(height=10); par(mfrow=c(3,3))
for(i in 1:m) {
  plot((fit_data$detected_v[,i]),xaxt='n',xlab='',ylab='new detected cases',main=paste0('vaccinated (',parms$ag_labels[i],')'),cex.main=1.25,cex.axis=1.25,cex.lab=1.25,
       ylim=c(0,max(model_out$detected_v_high))); 
  polygon(c(1:Tmax,Tmax:1),c(model_out$detected_v_high[,i],rev(model_out$detected_v_low[,i])),col='lightgray')
  points(fit_data$detected_v[,i])
  lines((model_out$detected_v[fit_period,i]),col='red',lwd=3)
  axis(1, at=seq(1,nrow(fit_data$detected_v),by=7), labels=format(dates[seq(1,nrow(fit_data$detected_v),by=7)],'%d/%m'), las=2 ,cex.axis=1.25)
}
x11(height=10); par(mfrow=c(3,3))
for(i in 1:m) {
  plot((fit_data$detected_b[,i]),xaxt='n',xlab='',ylab='new detected cases',main=paste0('boostered (',parms$ag_labels[i],')'),cex.main=1.25,cex.axis=1.25,cex.lab=1.25,
       ylim=c(0,max(model_out$detected_b_high))); 
  polygon(c(1:Tmax,Tmax:1),c(model_out$detected_b_high[,i],rev(model_out$detected_b_low[,i])),col='lightgray')
  points(fit_data$detected_b[,i])
  lines((model_out$detected_b[fit_period,i]),col='red',lwd=3)
  axis(1, at=seq(1,nrow(fit_data$detected_b),by=7), labels=format(dates[seq(1,nrow(fit_data$detected_b),by=7)],'%d/%m'), las=2 ,cex.axis=1.25)
}

x11(height=10); par(mfrow=c(3,3))
for(i in 1:m) {
  plot(fit_data$severe_nv[,i],xaxt='n',xlab='',ylab='new severe cases',main=paste0('unvaccinated (',parms$ag_labels[i],')'),cex.main=1.25,cex.axis=1.25,cex.lab=1.25,
       ylim=c(0,max(model_out$sevp_nv_high))); 
  polygon(c(1:Tmax,Tmax:1),c(model_out$sevp_nv_high[,i],rev(model_out$sevp_nv_low[,i])),col='lightgray')
  points(fit_data$severe_nv[,i])
  lines(model_out$sevp_nv[fit_period,i],col='red',lwd=3)
  axis(1, at=seq(1,nrow(fit_data$severe_nv),by=7), labels=format(dates[seq(1,nrow(fit_data$severe_nv),by=7)],'%d/%m'), las=2 ,cex.axis=1.25)
}
x11(height=10); par(mfrow=c(3,3))
for(i in 1:m) {
  plot((fit_data$severe_v[,i]),xaxt='n',xlab='',ylab='new severe cases',main=paste0('vaccinated (',parms$ag_labels[i],')'),cex.main=1.25,cex.axis=1.25,cex.lab=1.25,
       ylim=c(0,max(model_out$sevp_v_high))); 
  polygon(c(1:Tmax,Tmax:1),c(model_out$sevp_v_high[,i],rev(model_out$sevp_v_low[,i])),col='lightgray')
  points(fit_data$severe_v[,i])
  lines((model_out$sevp_v[fit_period,i]),col='red',lwd=3)
  axis(1, at=seq(1,nrow(fit_data$severe_v),by=7), labels=format(dates[seq(1,nrow(fit_data$severe_v),by=7)],'%d/%m'), las=2 ,cex.axis=1.25)
}
x11(height=10); par(mfrow=c(3,3))
for(i in 1:m) {
  plot((fit_data$severe_b[,i]),xaxt='n',xlab='',ylab='new severe cases',main=paste0('boostered (',parms$ag_labels[i],')'),cex.main=1.25,cex.axis=1.25,cex.lab=1.25,
       ylim=c(0,max(model_out$sevp_b_high))); 
  polygon(c(1:Tmax,Tmax:1),c(model_out$sevp_b_high[,i],rev(model_out$sevp_b_low[,i])),col='lightgray')
  points(fit_data$severe_b[,i])
  lines((model_out$sevp_b[fit_period,i]),col='red',lwd=3)
  axis(1, at=seq(1,nrow(fit_data$severe_b),by=7), labels=format(dates[seq(1,nrow(fit_data$severe_b),by=7)],'%d/%m'), las=2 ,cex.axis=1.25)
}

# Daily detected and severe per vac status - summary over all age-groups
x11(); par(mfrow=c(2,3))
plot(rowSums(fit_data$detected_nv),xaxt='n',xlab='',ylab='',main='detected_nv',cex.main=1.5,cex.axis=1.5,
     ylim=c(0,max(rowSums(model_out$detected_nv_high))))
polygon(c(fit_period,rev(fit_period)),c(rowSums(model_out$detected_nv_high[fit_period,]),rev(rowSums(model_out$detected_nv_low[fit_period,]))),col='lightgray')
points(rowSums(fit_data$detected_nv))
lines(smooth_window(rowSums(fit_data$detected_nv),7),lwd=3)
lines(rowSums(model_out$detected_nv[fit_period,]),col='red',lwd=3)
axis(1, at=seq(1,nrow(fit_data$detected_nv),by=7), labels=format(dates[seq(1,nrow(fit_data$detected_nv),by=7)],'%d/%m'), las=2 ,cex.axis=1.5)

plot(rowSums(fit_data$detected_v),xaxt='n',xlab='',ylab='',main='detected_v',cex.main=1.5,cex.axis=1.5,
     ylim=c(0,1000+max(rowSums(model_out$detected_v_high))))
polygon(c(fit_period,rev(fit_period)),c(rowSums(model_out$detected_v_high[fit_period,]),rev(rowSums(model_out$detected_v_low[fit_period,]))),col='lightgray')
points(rowSums(fit_data$detected_v))
lines(smooth_window(rowSums(fit_data$detected_v),7),lwd=3)
lines(rowSums(model_out$detected_v[fit_period,]),col='red',lwd=3)
axis(1, at=seq(1,nrow(fit_data$detected_v),by=7), labels=format(dates[seq(1,nrow(fit_data$detected_v),by=7)],'%d/%m'), las=2 ,cex.axis=1.5)

plot(rowSums(fit_data$detected_b),xaxt='n',xlab='',ylab='',main='detected_b',cex.main=1.5,cex.axis=1.5,
     ylim=c(0,1000+max(rowSums(model_out$detected_b_high))))
polygon(c(fit_period,rev(fit_period)),c(rowSums(model_out$detected_b_high[fit_period,]),rev(rowSums(model_out$detected_b_low[fit_period,]))),col='lightgray')
points(rowSums(fit_data$detected_b))
lines(smooth_window(rowSums(fit_data$detected_b),7),lwd=3)
lines(rowSums(model_out$detected_b[fit_period,]),col='red',lwd=3)
axis(1, at=seq(1,nrow(fit_data$detected_b),by=7), labels=format(dates[seq(1,nrow(fit_data$detected_b),by=7)],'%d/%m'), las=2 ,cex.axis=1.5)

plot(rowSums(fit_data$severe_nv),xaxt='n',xlab='',ylab='',main='severe_nv',cex.main=1.5,cex.axis=1.5,
     ylim=c(0,max(rowSums(model_out$sevp_nv_high))))
polygon(c(fit_period,rev(fit_period)),c(rowSums(model_out$sevp_nv_high[fit_period,]),rev(rowSums(model_out$sevp_nv_low[fit_period,]))),col='lightgray')
points(rowSums(fit_data$severe_nv))
lines(smooth_window(rowSums(fit_data$severe_nv),7),lwd=3)
lines(rowSums(model_out$sevp_nv[fit_period,]),col='red',lwd=3)
axis(1, at=seq(1,nrow(fit_data$severe_nv),by=7), labels=format(dates[seq(1,nrow(fit_data$severe_nv),by=7)],'%d/%m'), las=2 ,cex.axis=1.5)

plot(rowSums(fit_data$severe_v),xaxt='n',xlab='',ylab='',main='severe_v',cex.main=1.5,cex.axis=1.5,
     ylim=c(0,max(rowSums(model_out$sevp_v_high))))
polygon(c(fit_period,rev(fit_period)),c(rowSums(model_out$sevp_v_high[fit_period,]),rev(rowSums(model_out$sevp_v_low[fit_period,]))),col='lightgray')
points(rowSums(fit_data$severe_v))
lines(smooth_window(rowSums(fit_data$severe_v),7),lwd=3)
lines(rowSums(model_out$sevp_v[fit_period,]),col='red',lwd=3)
axis(1, at=seq(1,nrow(fit_data$severe_v),by=7), labels=format(dates[seq(1,nrow(fit_data$severe_v),by=7)],'%d/%m'), las=2 ,cex.axis=1.5)

plot(rowSums(fit_data$severe_b),xaxt='n',xlab='',ylab='',main='severe_b',cex.main=1.5,cex.axis=1.5,
     ylim=c(0,max(rowSums(model_out$sevp_b_high))))
polygon(c(fit_period,rev(fit_period)),c(rowSums(model_out$sevp_b_high[fit_period,]),rev(rowSums(model_out$sevp_b_low[fit_period,]))),col='lightgray')
points(rowSums(fit_data$severe_b))
lines(smooth_window(rowSums(fit_data$severe_b),7),lwd=3)
lines(rowSums(model_out$sevp_b[fit_period,]),col='red',lwd=3)
axis(1, at=seq(1,nrow(fit_data$severe_b),by=7), labels=format(dates[seq(1,nrow(fit_data$severe_b),by=7)],'%d/%m'), las=2 ,cex.axis=1.5)


x11(); barplot(parms$det_rates,ylab='reporting rates',names.arg=parms$ag_labels)


