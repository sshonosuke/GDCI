rm(list=ls())
set.seed(1)


## settings 
sc <- 2    # scenario number (from 1 to 4)
n <- 200    # number of samples
V <- 0.5     # random effect variance 
D <- rep(1, n)    # sampling variances 


## true signals 
if(sc==1){
  th <- sqrt(V)*rnorm(n)
}

if(sc==2){
  th <- sqrt(V)*rnorm(n)
  sub <- sample(1:n, n/20)
  th[sub] <- th[sub] + 7
}

if(sc==3){    
  th <- rlnorm(n, meanlog=0, sdlog=sqrt(V))
}

if(sc==4){    
  th <- rnorm(n, 0, 0.2)
  sub <- sample(1:n, n/10)
  th[sub] <- th[sub] + 10*V
}


## data generation 
Y <- th + sqrt(D)*rnorm(n)


## Proposed robust 95% interval 
source("GDCI-function.R")
sel <- GDCI.select(Y, D)
fit.GD <- GDCI(Y=Y, D=D, gam=sel$gam)
CI.GD <- fit.GD$CI


## Parametric 95% interval 
fit.EB <- GDCI(Y=Y, D=D, gam=0)
CI.EB <- fit.EB$CI

## AKM (95% interval)
library(ebci)
fit.AKM <- ebci(formula = Y~1, se = sqrt(D), weight=1/sqrt(D), alpha = 0.05)
CI.AKM <- cbind(fit.AKM$df$th_eb - fit.AKM$df$len_eb,fit.AKM$df$th_eb + fit.AKM$df$len_eb) 



## Result (coverage probability)
mean(CI.GD[,1]<th & CI.GD[,2]>th)
mean(CI.EB[,1]<th & CI.EB[,2]>th)
mean(CI.AKM[,1]<th & CI.AKM[,2]>th)


## Result (interval length)
mean(apply(CI.GD, 1, diff))
mean(apply(CI.EB, 1, diff))
mean(apply(CI.AKM, 1, diff))


