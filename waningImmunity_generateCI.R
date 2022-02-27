# model_out$detected_nv_low <- model_out$detected_nv - sqrt(model_out$detected_nv+model_out$detected_nv^2/parms$dp)
# model_out$detected_nv_high <- model_out$detected_nv + sqrt(model_out$detected_nv+model_out$detected_nv^2/parms$dp)
# model_out$detected_v_low <- model_out$detected_v - sqrt(model_out$detected_v+model_out$detected_v^2/parms$dp)
# model_out$detected_v_high <- model_out$detected_v + sqrt(model_out$detected_v+model_out$detected_v^2/parms$dp)
# model_out$detected_b_low <- model_out$detected_b - sqrt(model_out$detected_b+model_out$detected_b^2/parms$dp)
# model_out$detected_b_high <- model_out$detected_b + sqrt(model_out$detected_b+model_out$detected_b^2/parms$dp)
# 
# model_out$sevp_nv_low <- model_out$sevp_nv - sqrt(model_out$sevp_nv+model_out$sevp_nv^2/parms$dp)
# model_out$sevp_nv_high <- model_out$sevp_nv + sqrt(model_out$sevp_nv+model_out$sevp_nv^2/parms$dp)
# model_out$sevp_v_low <- model_out$sevp_v - sqrt(model_out$sevp_v+model_out$sevp_v^2/parms$dp)
# model_out$sevp_v_high <- model_out$sevp_v + sqrt(model_out$sevp_v+model_out$sevp_v^2/parms$dp)
# model_out$sevp_b_low <- model_out$sevp_b - sqrt(model_out$sevp_b+model_out$sevp_b^2/parms$dp)
# model_out$sevp_b_high <- model_out$sevp_b + sqrt(model_out$sevp_b+model_out$sevp_b^2/parms$dp)


getNB95pCI <- function(dat,dp) {
  n <- nrow(dat)
  m <- ncol(dat)
  ci95_low <- matrix(0,nrow=n,ncol=m)
  ci95_high <- matrix(0,nrow=n,ncol=m)
  for(i in 1:m) {
    obs <- sapply(1:1000,function(j) rnbinom(n=n,mu=dat[,i],size=dp))
    obss <- t(sapply(1:n, function(j) sort(obs[j,])))
    ci95_low[,i] <- obss[,51]
    ci95_high[,i] <- obss[,950]
  }
  return (list(low=ci95_low,high=ci95_high))
}

getPoisson95pCI <- function(dat) {
  n <- nrow(dat)
  m <- ncol(dat)
  ci95_low <- matrix(0,nrow=n,ncol=m)
  ci95_high <- matrix(0,nrow=n,ncol=m)
  for(i in 1:m) {
    obs <- sapply(1:1000,function(j) rpois(n=n,lambda=dat[,i]))
    obss <- t(sapply(1:n, function(j) sort(obs[j,])))
    ci95_low[,i] <- obss[,51]
    ci95_high[,i] <- obss[,950]
  }
  return (list(low=ci95_low,high=ci95_high))
}

if(fitUsingNB) {
  ci95_detected_nv <- getNB95pCI(model_out$detected_nv,parms$dp)
  ci95_detected_v  <- getNB95pCI(model_out$detected_v,parms$dp)
  ci95_detected_b  <- getNB95pCI(model_out$detected_b,parms$dp)
  ci95_sevp_nv     <- getNB95pCI(model_out$sevp_nv,parms$dp)
  ci95_sevp_v      <- getNB95pCI(model_out$sevp_v,parms$dp)
  ci95_sevp_b      <- getNB95pCI(model_out$sevp_b,parms$dp)
} else {
  ci95_detected_nv <- getPoisson95pCI(model_out$detected_nv)
  ci95_detected_v  <- getPoisson95pCI(model_out$detected_v)
  ci95_detected_b  <- getPoisson95pCI(model_out$detected_b)
  ci95_sevp_nv     <- getPoisson95pCI(model_out$sevp_nv)
  ci95_sevp_v      <- getPoisson95pCI(model_out$sevp_v)
  ci95_sevp_b      <- getPoisson95pCI(model_out$sevp_b)
}



model_out$detected_nv_low <- ci95_detected_nv$low
model_out$detected_nv_high <- ci95_detected_nv$high
model_out$detected_v_low <- ci95_detected_v$low
model_out$detected_v_high <- ci95_detected_v$high
model_out$detected_b_low <- ci95_detected_b$low
model_out$detected_b_high <- ci95_detected_b$high

model_out$sevp_nv_low  <- ci95_sevp_nv$low
model_out$sevp_nv_high <- ci95_sevp_nv$high
model_out$sevp_v_low   <- ci95_sevp_v$low
model_out$sevp_v_high  <- ci95_sevp_v$high
model_out$sevp_b_low   <- ci95_sevp_b$low
model_out$sevp_b_high  <- ci95_sevp_b$high