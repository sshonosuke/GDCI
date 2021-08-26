###  robust empirical Bayes confidence intervals using gamma-divergence (given gamma)
# Y: response vector 
# D: vector of sampling variance 
# X: design matrix (default is NULL)
# gam: tuning parameter 
# alpha: significance level (default is 0.05) 
# ad: adjustment factor (default is 0)

GDCI <- function(Y, D, X=NULL, gam=0.1, alpha=0.05, ad=0){
  # preparation
  n <- length(Y)
  if(is.null(X)){ X <- matrix(1, n, 1) }
  p <- dim(X)[2]
  bb <- 10^10     # upper bound of parameters 
  ep <- 10^(-10)    # lower bound of variance parameter 
  zz <- qnorm(1-alpha/2)
  
  # initial value (standard parametric method) 
  Q <- function(para){
    mu <- as.vector(X%*%para[1:p])
    A <- para[p+1]
    val <- sum(log(A+D)) + sum((Y-mu)^2/(A+D)) - ad*log(A) 
    if(is.na(val)){ val <- 10^10 }
    return(val)
  }
  ML <- optim(fn=Q, par=c(rep(0, p), 1), method="L-BFGS-B", lower=c(rep(-bb, p), ep), upper=rep(bb, p+1), hessian=T)
  hbeta <- ML$par[1:p]
  hA <- ML$par[p+1]
  
  # fitting 
  if(gam==0){
    hmu <- as.vector(X%*%hbeta)
    EB <- Y - D/(hA+D)*(Y-hmu) 
    ss <- sqrt( hA*D/(hA+D) )
    CI <- cbind(EB-zz*ss, EB+zz*ss)
    Result <- list(CI=CI, EB=EB, mu=hmu, beta=hbeta, tau2=hA, posvar=ss^2)
  }
  if(gam>0){
    Q <- function(para){
      mu <- as.vector(X%*%para[1:p])
      A <- para[p+1]
      cc <- (2*pi*(A+D))^(0.5*gam^2/(1+gam)) 
      GD <- cc*dnorm(Y, mu, sqrt(A+D))^gam/gam 
      val <- sum(-GD) - ad*log(A)
      if(is.na(val)){ val <- 10^10 }
      return(val)
    }
    RML <- optim(fn=Q, par=c(hbeta, hA), method="L-BFGS-B", lower=c(rep(-bb, p), ep), upper=rep(bb, p+1), hessian=T)
    rbeta <- RML$par[1:p]
    rmu <- as.vector(X%*%rbeta)
    rA <- RML$par[p+1]
    ww <- dnorm(Y, rmu, sqrt(rA+D))^gam
    cc <- (2*pi*(rA+D))^(0.5*gam^2/(1+gam)) 
    rEB <- Y - D/(rA+D)*(Y-rmu)*ww*cc 
    rss2 <- D + cc*ww*D^2/(rA+D)^2*(gam*(Y-rmu)^2 - (rA+D))
    rss2[rss2<ep] <- ep
    rss <-  sqrt(rss2)
    rCI <- cbind(rEB-zz*rss, rEB+zz*rss)
    Result <- list(CI=rCI, EB=rEB, mu=rmu, beta=rbeta, tau2=rA, posvar=rss2)
  }
  return(Result)
}








###  tuning parameter selection for GDCI
# Y: response vector 
# D: vector of sampling variance 
# X: design matrix (default is NULL)
# gam.set: set of candidate values of tuning parameter 

GDCI.select <- function(Y, D, X=NULL, gam.set=NULL, ad=0){
  if(is.null(gam.set)){  gam.set <- seq(0, 0.5, by=0.01)  }
  L <- length(gam.set)
  PV <- c()
  for(k in 1:L){
    PV[k] <- mean( GDCI(Y=Y, D=D, X=X, gam=gam.set[k], ad=ad)$posvar/D )
  }
  opt.gam <- gam.set[which.min(PV)]
  Result <- list(gam=opt.gam, PV=PV)
  return(Result)
}









