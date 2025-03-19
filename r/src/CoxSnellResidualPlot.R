# functions for plotting Cox-Snell Residuals

plotCoxSnellCPH<- function(survFit, delta, fitType){
  
  # get Cox-Snell Residual based on Martingale Residuals
  mg.residual <- resid(survFit, type = "martingale")
  
  cs.residual <- delta - mg.residual
  
  # Graphical Plot 
  fit.cs <- survfit(Surv(cs.residual, delta) ~ 1) # get KM estimates
  H.cs <- cumsum(fit.cs$n.event/fit.cs$n.risk)
  
  plot(fit.cs$time, H.cs, type='s', col='blue', 
       main = paste0('Cox-PH - ', fitType), 
       xlab = 'Residual', ylab = 'Nelson-Aalen Cum. Hazard') 
  #Note here that 'time' is the value of the Cox-Snell residual
  abline(0, 1, col='red',  lty = 2)
}

plotCSExpAFT <- function(survFit, data, fitType){
  sigma <- survFit$scale
  eta   <- -survFit$linear.predictors/sigma
  
  r.exp <- data$time * exp(eta)
  
  fit   <- survfit(Surv(r.exp, data$delta) ~ 1)
  H.exp <- cumsum(fit$n.event / fit$n.risk)
  
  plot(H.exp ~ fit$time, type = 'l', main = paste0('Exponential AFT - ', fitType),
       ylab = 'Estimated Cumulative Hazard', xlab = 'Cox-Snell Residual')
  abline(0, 1, col='red',  lty=2)
}

plotCSWeibullAFT <- function(survFit, data, fitType){
  sigma  <- survFit$scale
  alpha  <- 1 / sigma
  eta    <- -survFit$linear.predictors / sigma
  
  r.wb <- data$time^alpha * exp(eta)
  
  fit   <- survfit(Surv(r.wb, data$delta) ~ 1)
  H.wb  <- cumsum(fit$n.event/fit$n.risk)
  
  plot(H.wb ~ fit$time, type = 'l', main = paste0('Weibull AFT - ', fitType),
       ylab = 'Estimated Cumulative Hazard', xlab = 'Cox-Snell Residual')
  abline(0, 1, col='red',  lty=2)
}

plotCSLogLogisticAFT <- function(survFit, data, fitType){
  sigma  <- survFit$scale
  alpha  <- 1 / sigma
  eta    <- -survFit$linear.predictors / sigma
  
  r.ll  <- -log(1/(1 + data$time^alpha*exp(eta)))
  
  fit   <- survfit(Surv(r.ll, data$delta) ~ 1)
  H.ll  <- cumsum(fit$n.event / fit$n.risk)
  
  plot(H.ll ~ fit$time, type = 'l', main = paste0('Log-Logistic AFT- ', fitType),
       ylab = 'Estimated Cumulative Hazard', xlab = 'Cox-Snell Residual')
  abline(0, 1, col='red',  lty=2)
}

plotCSLogNormalAFT <- function(survFit, data, fitType){
  eta    <- -survFit$linear.predictors / survFit$scale
  r.ln   <- -log(1 - pnorm((log(data$time) - survFit$linear.predictors) / survFit$scale))
  
  fit   <- survfit(Surv(r.ln, data$delta) ~ 1)
  H.ln  <- cumsum(fit$n.event/fit$n.risk)
  
  plot(H.ln ~ fit$time, type = 'l', main = paste0('Log-Normal AFT - ', fitType),
       ylab = 'Estimated Cumulative Hazard', xlab = 'Cox-Snell Residual')
  abline(0, 1, col='red',  lty = 2) 
}
