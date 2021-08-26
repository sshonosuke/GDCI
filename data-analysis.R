rm(list=ls())

## load dataset 
load("crimedata.RData")

## The dataset includes the following three objects. 
# mY: vector of average counts 
# X: matrix of covariates (including an intercept term)
# SV: vector of sampling variances



## Proposed robust 95% interval 
source("GDCI-function.R")

gam.set <- seq(0, 0.2, by=0.005)
sel <- GDCI.select(Y=mY, D=SV, X=X, gam.set=gam.set)
sel$gam
fit.GD <- GDCI(Y=mY, D=SV, X=X, gam=sel$gam)
fit.GD$beta     # regression coefficients
CI.GD <- fit.GD$CI


## Parametric 95% intervals
fit.EB <- GDCI(Y=mY, D=SV, X=X, gam=0)
fit.EB$beta     # regression coefficients
CI.EB <- fit.EB$CI


## AKM method (95% interval)
library(ebci)
fit.AKM <- ebci(formula=mY~X[,-1], se=sqrt(SV), weight=1/sqrt(SV), alpha=0.05)
fit.AKM$delta    # regression coefficients
CI.AKM <- cbind(fit.AKM$df$th_eb - fit.AKM$df$len_eb, fit.AKM$df$th_eb + fit.AKM$df$len_eb) 



## ML 95% interval (direct intervals using observations and their variances)
CI.ML <- cbind(mY+qnorm(0.025)*sqrt(SV), mY+qnorm(0.975)*sqrt(SV))



## interval length
len1 <- apply(CI.GD, 1, diff)
len2 <- apply(CI.EB, 1, diff)
len3 <- apply(CI.AKM, 1, diff)
len4 <- apply(CI.ML, 1, diff)




###   Figure 1 in the manuscript   ###
par(mfcol=c(1,2))
# selection criterion as a function of gamma
plot(gam.set, sel$PV, type="l", ylab="criterion", xlab="gamma")
# interval length
ran <- range(len1/len4, len2/len4, len3/len4)
sY <- mY/sqrt(SV)
plot(sY, len1/len4, log="x", ylim=ran, pch=8, cex=0.8, 
     ylab="Ratio of interval length (to ML)", xlab="Scaled observed value")
points(sY, len2/len4, col="red", pch=1, cex=0.8)
points(sY, len3/len4, col="blue", pch=16, cex=0.8)
legend("bottomleft", legend=c("GD (proposed)", "EB", "AKM"), 
       col=c("black", "red", "blue"), pch=c(8,1,16))




###   Figure 2 in the manuscript    ###
par(mfcol=c(1,2))
# large signals 
S <- 30
sub <- order(mY, decreasing=T)[1:S]
mat <- cbind(mY[sub], CI.GD[sub,], CI.EB[sub,], CI.AKM[sub,], CI.ML[sub,])
ran <- range(mat)
bb <- c(-0.1, -0.1, 0.1, 0.1)
cc <- 0.1
plot((1:S)+bb[4], mat[,1], pch=8, ylim=ran, col="red", ylab="crime number", 
     xlab="area", main="outlying observations")
for(k in 1:S){
  lines(rep(k+bb[1],2), mat[k,2:3], col=1)
  lines(bb[1]+c(k-cc, k+cc), rep(mat[k,2],2), col=1)
  lines(bb[1]+c(k-cc, k+cc), rep(mat[k,3],2), col=1)
  lines(rep(k+bb[2],2), mat[k,4:5], col="blue")
  lines(bb[2]+c(k-cc, k+cc), rep(mat[k,4],2), col="blue")
  lines(bb[2]+c(k-cc, k+cc), rep(mat[k,5],2), col="blue")
  lines(rep(k+bb[3],2), mat[k,6:7], col="blue", lty=2)
  lines(bb[3]+c(k-cc, k+cc), rep(mat[k,6],2), col="blue", lty=2)
  lines(bb[3]+c(k-cc, k+cc), rep(mat[k,7],2), col="blue", lty=2)
  lines(rep(k+bb[4],2), mat[k,8:9], col="red")
  lines(bb[4]+c(k-cc, k+cc), rep(mat[k,8],2), col="red")
  lines(bb[4]+c(k-cc, k+cc), rep(mat[k,9],2), col="red")
}
legend("topright", legend=c("GD", "EB", "AKM", "ML"), lty=c(1,1,2,1), 
       col=c("black", "blue", "blue", "red"), pch=c(NA, NA, NA, 8))

# moderate signals (figure)
sub <- which(mY/sqrt(SV)<5 & mY/sqrt(SV)>1.5 & (len1/len3)<0.9)
S <- 30
sub <- sample(sub, S)
sub <- sub[order(mY[sub], decreasing=T)]
mat <- cbind(mY[sub], CI.GD[sub,], CI.ML[sub,])
ran <- range(mat)
bb <- c(-0.2, 0)
cc <- 0.1
plot(1:S, mat[,1], pch=8, xlab="area", ylim=ran, ylab="crime number", 
     col="red", main="non-outlying observations")
for(k in 1:S){
  lines(rep(k+bb[1],2), mat[k,2:3], col=1)
  lines(bb[1]+c(k-cc, k+cc), rep(mat[k,2],2), col=1)
  lines(bb[1]+c(k-cc, k+cc), rep(mat[k,3],2), col=1)
  lines(rep(k+bb[2],2), mat[k,4:5], col="red")
  lines(bb[2]+c(k-cc, k+cc), rep(mat[k,4],2), col="red")
  lines(bb[2]+c(k-cc, k+cc), rep(mat[k,5],2), col="red")
}
legend("topright", legend=c("GD", "ML"), lty=1, col=c("black", "red"), pch=c(NA, 8))




