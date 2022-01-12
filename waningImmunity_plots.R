
x11(width=10); matplot(mobility_data[fit_period,],type='l',lty=1:4,lwd=3,xaxt='n',cex.axis=1.5,ylab='1-(% change from baseline)/100 ',cex.axis=1.5,cex.lab=1.5)
axis(1, at=seq(1,fit_length,by=14), labels=format(dates[seq(1,fit_length,by=14)],'%d/%m'), las=2 ,cex.axis=1.5)
legend('bottomright',c('household','work','school','community'),lwd=rep(3,4),lty=c(1,2,4,3),col=c(1,2,4,3),cex=1.5)

x11(); plot(model_out$Reff,col='blue',lwd=3,type='l',xaxt='n',xlab='',ylab='R',ylim=c(0.7,1.5),cex.axis=1.5,cex.lab=1.5);
# lines(model_out$Reff_observed)
axis(1, at=seq(1,length(model_out$Reff_observed),by=7), labels=format(dates[seq(1,length(model_out$Reff_observed),by=7)],'%d/%m'), las=2 ,cex.axis=1.5)
grid()


# Daily detected and severe per vac status and age-group

x11(height=10); par(mfrow=c(3,3))
for(i in 1:m) {
  plot(fit_data$detected_nv[,i],xaxt='n',xlab='',ylab='new detected cases',main=paste0('unvaccinated (',parms$ag_labels[i],')'),cex.main=1.25,cex.axis=1.25,cex.lab=1.25,
       ylim=c(0,max(fit_data$detected_nv))); 
  # lines(model_out$detected_nv[,i],col='blue',lwd=3)
  lines(model_out$detected_nv[fit_period,i],col='red',lwd=3)
  axis(1, at=seq(1,nrow(fit_data$detected_nv),by=7), labels=format(dates[seq(1,nrow(fit_data$detected_nv),by=7)],'%d/%m'), las=2 ,cex.axis=1.25)
}
x11(height=10); par(mfrow=c(3,3))
for(i in 1:m) {
  plot((fit_data$detected_v[,i]),xaxt='n',xlab='',ylab='new detected cases',main=paste0('vaccinated (',parms$ag_labels[i],')'),cex.main=1.25,cex.axis=1.25,cex.lab=1.25,
       ylim=c(0,max(fit_data$detected_v))); 
  # lines((model_out$detected_v[,i]),col='blue',lwd=3)
  lines((model_out$detected_v[fit_period,i]),col='red',lwd=3)
  axis(1, at=seq(1,nrow(fit_data$detected_v),by=7), labels=format(dates[seq(1,nrow(fit_data$detected_v),by=7)],'%d/%m'), las=2 ,cex.axis=1.25)
}
x11(height=10); par(mfrow=c(3,3))
for(i in 1:m) {
  plot((fit_data$detected_b[,i]),xaxt='n',xlab='',ylab='new detected cases',main=paste0('boostered (',parms$ag_labels[i],')'),cex.main=1.25,cex.axis=1.25,cex.lab=1.25,
       ylim=c(0,max(fit_data$detected_b))); 
  # lines((model_out$detected_b[,i]),col='blue',lwd=3)
  lines((model_out$detected_b[fit_period,i]),col='red',lwd=3)
  axis(1, at=seq(1,nrow(fit_data$detected_b),by=7), labels=format(dates[seq(1,nrow(fit_data$detected_b),by=7)],'%d/%m'), las=2 ,cex.axis=1.25)
}

x11(height=10); par(mfrow=c(3,3))
for(i in 1:m) {
  plot(fit_data$severe_nv[,i],xaxt='n',xlab='',ylab='new severe cases',main=paste0('unvaccinated (',parms$ag_labels[i],')'),cex.main=1.25,cex.axis=1.25,cex.lab=1.25,
       ylim=c(0,max(fit_data$severe_nv))); 
  # lines(model_out$sevp_nv[,i],col='blue',lwd=3)
  lines(model_out$sevp_nv[fit_period,i],col='red',lwd=3)
  axis(1, at=seq(1,nrow(fit_data$severe_nv),by=7), labels=format(dates[seq(1,nrow(fit_data$severe_nv),by=7)],'%d/%m'), las=2 ,cex.axis=1.25)
}
x11(height=10); par(mfrow=c(3,3))
for(i in 1:m) {
  plot((fit_data$severe_v[,i]),xaxt='n',xlab='',ylab='new severe cases',main=paste0('vaccinated (',parms$ag_labels[i],')'),cex.main=1.25,cex.axis=1.25,cex.lab=1.25,
       ylim=c(0,max(fit_data$severe_nv))); 
  # lines((model_out$sevp_v[,i]),col='blue',lwd=3)
  lines((model_out$sevp_v[fit_period,i]),col='red',lwd=3)
  axis(1, at=seq(1,nrow(fit_data$severe_v),by=7), labels=format(dates[seq(1,nrow(fit_data$severe_v),by=7)],'%d/%m'), las=2 ,cex.axis=1.25)
}
x11(height=10); par(mfrow=c(3,3))
for(i in 1:m) {
  plot((fit_data$severe_b[,i]),xaxt='n',xlab='',ylab='new severe cases',main=paste0('boostered (',parms$ag_labels[i],')'),cex.main=1.25,cex.axis=1.25,cex.lab=1.25,
       ylim=c(0,max(fit_data$severe_nv))); 
  # lines((model_out$sevp_b[,i]),col='blue',lwd=3)
  lines((model_out$sevp_b[fit_period,i]),col='red',lwd=3)
  axis(1, at=seq(1,nrow(fit_data$severe_b),by=7), labels=format(dates[seq(1,nrow(fit_data$severe_b),by=7)],'%d/%m'), las=2 ,cex.axis=1.25)
}

# Daily detected and severe per vac status - summary over all age-groups
x11(); par(mfrow=c(2,3))
plot(rowSums(fit_data$detected_nv),xaxt='n',xlab='',ylab='',main='detected_nv',cex.main=1.5,cex.axis=1.5,
     ylim=c(0,max(rowSums(fit_data$detected_nv[,]))))
lines(smooth_window(rowSums(fit_data$detected_nv),7),lwd=3)
lines(rowSums(model_out$detected_nv),col='blue',lwd=3)
lines(rowSums(model_out$detected_nv[fit_period,]),col='red',lwd=3)
axis(1, at=seq(1,nrow(fit_data$detected_nv),by=7), labels=format(dates[seq(1,nrow(fit_data$detected_nv),by=7)],'%d/%m'), las=2 ,cex.axis=1.5)

plot(rowSums(fit_data$detected_v),xaxt='n',xlab='',ylab='',main='detected_v',cex.main=1.5,cex.axis=1.5,
     ylim=c(0,1000+max(rowSums(fit_data$detected_v[,]))))
lines(smooth_window(rowSums(fit_data$detected_v),7),lwd=3)
lines(rowSums(model_out$detected_v),col='blue',lwd=3)
lines(rowSums(model_out$detected_v[fit_period,]),col='red',lwd=3)
axis(1, at=seq(1,nrow(fit_data$detected_v),by=7), labels=format(dates[seq(1,nrow(fit_data$detected_v),by=7)],'%d/%m'), las=2 ,cex.axis=1.5)

plot(rowSums(fit_data$detected_b),xaxt='n',xlab='',ylab='',main='detected_b',cex.main=1.5,cex.axis=1.5,
     ylim=c(0,1000+max(rowSums(fit_data$detected_b[,]))))
lines(smooth_window(rowSums(fit_data$detected_b),7),lwd=3)
lines(rowSums(model_out$detected_b),col='blue',lwd=3)
lines(rowSums(model_out$detected_b[fit_period,]),col='red',lwd=3)
axis(1, at=seq(1,nrow(fit_data$detected_b),by=7), labels=format(dates[seq(1,nrow(fit_data$detected_b),by=7)],'%d/%m'), las=2 ,cex.axis=1.5)

plot(rowSums(fit_data$severe_nv[,]),xaxt='n',xlab='',ylab='',main='severe_nv',cex.main=1.5,cex.axis=1.5,
     ylim=c(0,max(rowSums(fit_data$severe_nv))))
lines(smooth_window(rowSums(fit_data$severe_nv),7),lwd=3)
lines(rowSums(model_out$sevp_nv),col='blue',lwd=3)
lines(rowSums(model_out$sevp_nv[fit_period,]),col='red',lwd=3)
axis(1, at=seq(1,nrow(fit_data$severe_nv),by=7), labels=format(dates[seq(1,nrow(fit_data$severe_nv),by=7)],'%d/%m'), las=2 ,cex.axis=1.5)

plot(rowSums(fit_data$severe_v[,]),xaxt='n',xlab='',ylab='',main='severe_v',cex.main=1.5,cex.axis=1.5,
     ylim=c(0,max(rowSums(fit_data$severe_v))))
lines(smooth_window(rowSums(fit_data$severe_v),7),lwd=3)
lines(rowSums(model_out$sevp_v),col='blue',lwd=3)
lines(rowSums(model_out$sevp_v[fit_period,]),col='red',lwd=3)
axis(1, at=seq(1,nrow(fit_data$severe_v),by=7), labels=format(dates[seq(1,nrow(fit_data$severe_v),by=7)],'%d/%m'), las=2 ,cex.axis=1.5)

plot(rowSums(fit_data$severe_b[,]),xaxt='n',xlab='',ylab='',main='severe_b',cex.main=1.5,cex.axis=1.5,
     ylim=c(0,max(rowSums(fit_data$severe_b))))
lines(smooth_window(rowSums(fit_data$severe_b),7),lwd=3)
lines(rowSums(model_out$sevp_b),col='blue',lwd=3)
lines(rowSums(model_out$sevp_b[fit_period,]),col='red',lwd=3)
axis(1, at=seq(1,nrow(fit_data$severe_b),by=7), labels=format(dates[seq(1,nrow(fit_data$severe_b),by=7)],'%d/%m'), las=2 ,cex.axis=1.5)


# Daily detected - total
x11()
plot(rowSums(model_out$detected),main='Daily detected',type='l',col='red',lwd=3,xlab='',ylab='',xaxt='n',cex.axis=1.5,ylim=c(0,12000))
points(rowSums(fit_data$detected_nv+fit_data$detected_v+fit_data$detected_b))
axis(1, at=seq(1,fit_length,by=7), labels=format(fit_dates[seq(1,fit_length,by=7)],'%d/%m'), las=2 ,cex.axis=1.5)
legend('topleft',legend=c('data','simulation'),
       col=c('black','red'),lty=c(NA,1),pch=c(1,NA),lwd=c(NA,3),cex=1.5)


# Daily severe - total
x11()
plot(rowSums(model_out$sevp),main='Daily severe',type='l',col='red',lwd=3,xlab='',ylab='',xaxt='n',cex.axis=1.5,ylim=c(0,120))
points(rowSums(fit_data$severe_nv+fit_data$severe_v+fit_data$severe_b))
axis(1, at=seq(1,fit_length,by=7), labels=format(fit_dates[seq(1,fit_length,by=7)],'%d/%m'), las=2 ,cex.axis=1.5)
legend('topleft',legend=c('data','simulation'),
       col=c('black','red'),lty=c(NA,1),pch=c(1,NA),lwd=c(NA,3),cex=1.5)


x11(); barplot(1-(S0_nv+colSums(S0_v))/N,ylab='initial attack rates',names.arg=parms$ag_labels)


