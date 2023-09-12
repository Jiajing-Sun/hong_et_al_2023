################################################################################
##################################### FUNCTIONS ################################ 
################################################################################

# this function does the following: 
# 1. interpolate the dataset if it is at the beginning of the series 
# 2. setting the missing values at the beginning of the series to zero. 

clean.data.fcn <- function (x )
{
  x <- as.numeric(x) 
  temp0<-na_interpolation(x)
  temp<- which(is.na(x)==TRUE)
  temp1<- split(temp, cumsum(c(1, diff(temp) != 1)))
  if ((temp1$"1")[1]==1) { temp0 [temp1$"1" ] <-0 }
  x <- temp0
} 
################################################################################
# return series 

price2return<- function ( p )
{
  N<-length(p)
  return.series<- (log(p [2:N ]) -log(p[1:(N -1)]))*100 
  return.series<-na_interpolation(return.series)
  return (return.series )
}
################################################################################

specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k)) 
################################################################################
# function: generate lagged y 

lagit4me <- function(lag){ # getting lagged y, replace NA with zeros. 
  n <- length(y);
  pad <- rep(0,lag);
  return(c(pad,y)[1:n]);
} 
################################################################################
# normality test of Bai and Ng (2005) JBES 

bai.ng.normality <- function( x ) { 
  miu <- mean(x); sigma.2 <- var(x); miu.3 <-mean ((x-miu)^3); sigma.3 <- (mean((x-miu)^2))^(3/2)
  tau <- miu.3 / sigma.3# skewness
  alpha <- c(0, -3*sigma.2,-3* sqrt(sigma.2) * tau/2 )
  Z <- cbind((x-miu)^3 - miu.3* rep(1, length(x)), (x-miu* rep(1, length(x))), (x-miu)^2 - sigma.2 * rep(1, length(x))) 
  # capital.tau <- cov( Z)
  # use kernel based estimator for capital.tau as suggested on page 51 of Bai and Ng (2005)
  capital.tau <- lrvar (Z) * length(x) 
  capital.tau.22 <-capital.tau[1:2, 1:2] 
  alpha.2 <- c(1, -3* sigma.2) 
  s.miu.3 <- sqrt (t(alpha.2)%*% capital.tau.22%*% (alpha.2)) 
  pi.3 <- sqrt(length(x) )* miu.3 / (s.miu.3) 
  miu.4 <-mean ((x-miu)^4); sigma.4 <- (mean((x-miu)^2))^2 
  kappa <- miu.4 / sigma.4 
  beta <- c (1, -4*miu.3, -2*sigma.2*kappa ) 
  W<- cbind((x-miu)^4 - miu.4* rep(1, length(x)), (x-miu* rep(1, length(x))), (x-miu)^2 - sigma.2 * rep(1, length(x))) 
  capital.omega <- lrvar (W) * length(x) 
  s.kappa<- sqrt (t(beta)%*% capital.omega %*% (beta) / (sigma.2 ^ 4)) 
  pi.4 <- sqrt(length(x) )* (kappa-3) / s.kappa 
  bai.ng.2005<- pi.3^2 + pi.4^2 
  p.bai.ng.2005<- pchisq(bai.ng.2005, df=2, lower.tail = F)
  return( list (statistic=bai.ng.2005 ,p.value =p.bai.ng.2005))
}
################################################################################
# summary statistics 
summary.statistics <-function(x){ 
  x<-na.omit(x)
  min.x<- min(x)
  max.x<- max(x)
  average.x<- mean(x)
  sd.x <- sqrt(mean((x- average.x )^2)) 
  skewness.x<- skewness(x)
  kurtosis<- kurtosis(x)
  adf.stat<- adf.test(x)$statistic
  adf.p<- adf.test(x)$p.value
  bai.ng.2005.stat<-bai.ng.normality(x)$statistic
  bai.ng.2005.p<-bai.ng.normality(x)$p.value 
  result<- cbind(min.x, max.x, average.x, sd.x, skewness.x, kurtosis, 
                 adf.stat, adf.p , bai.ng.2005.stat,bai.ng.2005.p ) 
  return(result)
} 
################################################################################

s_t1_t2.fcn <- function( t1=t1 , t2=t2 ) 
{ if (t1<=t2) {
  temp <- sum(y [ t1:t2] ) 
} else {
  temp <-0
} 
  return(temp)
  
} 

################################################################################
# filtering standardized returns. 

standardized.r <- function( y = y ) { 
  library(rugarch) 
  spec = ugarchspec(mean.model = list(armaOrder = c(1,1), include.mean = TRUE),
                    variance.model = list(model = "sGARCH"), distribution.model = "norm")
  
  garch <- ugarchfit(spec = spec, data = y, drop = FALSE) 
  sigma.2 <- garch@fit$var 
  if (is.null (sigma.2)== TRUE) { standardized.r <- y / sqrt(var(y)) } else { 
    standardized.r <- y/sqrt(sigma.2 ) }
  return( standardized.r =standardized.r ) 
}
################################################################################
# arma(1,1) residuals 

arma.11.residuals <- function(y) { 
  temp <- arma(y , order = c(1, 1) )
  arma11.resid <- na.omit ( temp$residuals)
}


################################################################################
# CUSUM process for correlation coefficient. 

T.n.k.tilda.fcn <- function( W.vec ) { 
  sample.size.v2 <- length( W.vec ) 
  times <-( 1: sample.size.v2 ) 
  temp <- (W.vec - W.vec [sample.size.v2] )* times
}


T.n.k.tilda.fcn.v2 <- function( W.vec ) { 
  sample.size.v2 <- length( W.vec ) 
  times <-( 1: sample.size.v2 ) 
  temp <- (W.vec )* times
}

################################################################################
# generate correlation coefficients for G test statistic. 

corr.mat.fcn <-function(input.matrix= a.t , t1=1, t2= k ) { 
  temp.ind <- upper.tri ( diag ( dim(input.matrix)[2])) 
  temp1<- input.matrix[t1:t2, ] 
  if ( length( temp1) == dim(input.matrix)[2] ) { temp2 <- ( rep(1, dim(input.matrix)[2] * ( dim(input.matrix)[2] -1)/2 ))
  } else { temp2 <- cor(input.matrix[t1:t2,]) [ temp.ind] }
  return(temp2) 
}

################################################################################
# adjusted-range based KS test statistic 


ks.test.fcn <- function(res=res) { 
  res =na.omit(res)
  sample.size <- length(res)
  average.res <-mean (res ) 
  demeaned.res <-res - average.res 
  cumsum.range <- cumsum(demeaned.res)
  R.hat <- ( max(cumsum.range) - min (cumsum.range))/sqrt( sample.size) 
  KS.R <- (max(abs(cumsum.range))/sqrt( sample.size) )/R.hat # adjusted-range based KS test statistic 
  
  library(cointReg) 
  
  KS.numerator <-(max(abs(cumsum.range)) /sqrt(sample.size)) 
  s.e.Bartlett.denominator <- getLongRunVar(res, kernel = "ba", bandwidth = "and")
  s.e.Bartlett.denominator <- sqrt(s.e.Bartlett.denominator $Omega) 
  KS.0 <- KS.numerator / s.e.Bartlett.denominator # KS statistic based on small-b asymptotics 
  
  return (list (KS.R= KS.R, KS.0=KS.0 ) ) 
} 


################################################################################
# Gn test statistic of Shao and Zhang (2010)

Gn.univariate <- function( res = res ){ 
  s_t1_t2.fcn <- function( t1=t1 , t2=t2 ) { # function : forward and backward sum 
    if (t1<=t2) {
      temp <- sum(res [ t1:t2] ) 
    } else {
      temp <-0
    } 
  } 
  
  
  res =na.omit(res) 
  sample.size <- length(res)
  average.res <-mean (res ) 
  demeaned.res <-res - average.res 
  v.k <- rep( 0, sample.size-1) 
  Gn <- rep ( 0, sample.size-1) 
  for ( k in 1: (sample.size-1))
  { 
    temp1 <-0 
    for ( j in 1:k) { temp1 <- temp1+ (s_t1_t2.fcn(t1=1, t2=j ) -j/k* s_t1_t2.fcn(t1=1, t2=k ) ) ^2 }
    temp2 <-0 
    for ( j in (k+1): (sample.size-1) ) { temp2 <- temp2 + (s_t1_t2.fcn(t1=j, t2=( sample.size-1) ) - (sample.size - j+1)/(sample.size - k)* s_t1_t2.fcn(t1= (k+1), t2=(sample.size-1) )) ^2 } 
    v.k [k] <- (temp1+ temp2)/ sample.size^2 
    demeaned.res.v2 <- demeaned.res [ 1: (sample.size-1)]
    temp3 <- cumsum(demeaned.res.v2) /sqrt(sample.size-1)
    Gn[k] <- temp3 [k]^2 / v.k [k] 
  } 
  Gn <- max(Gn) 
  return(Gn=Gn)
}

################################################################################
# correlation matrix test for choi and shin (2020) 
choi.shin.cor.fcn <- function ( cor.vector =cor.vector ){ 
  s_t1_t2.fcn <- function( t1=t1 , t2=t2 ) 
  { if (t1<=t2) {
    temp <- sum(y [ t1:t2] ) 
  } else {
    temp <-0
  } 
    return(temp)
  } 
  
  y<- rowSums(abs( cor.vector) )
  sample.size <- length(y)
  
  Vk <- rep( 0, sample.size-1) 
  Gn.statistic <- rep ( 0, sample.size-1) 
  
  
  for ( k in 1: (sample.size-1))
  { 
    temp1 <-0 
    for ( j in 1:k) { temp1 <- temp1+ (s_t1_t2.fcn(t1=1, t2=j ) -j/k* s_t1_t2.fcn(t1=1, t2=k ) ) ^2 }
    temp2 <-0 
    for ( j in (k+1): (sample.size-1) ) { temp2 <- temp2 + (s_t1_t2.fcn(t1=j, t2= sample.size-1 ) - (sample.size - j+1)/(sample.size - k)* s_t1_t2.fcn(t1= (k+1), t2=sample.size-1 ) ) ^2 } 
    Vk [k] <- (temp1+ temp2)/ sample.size^2 
    
    #demeaned.y.v2 <- demeaned.y [ 1: (sample.size-1)] 
    
    demeaned.y.v2 <- y [ 1: (sample.size-1)] - y [sample.size]
    temp3 <- cumsum(demeaned.y.v2) /sqrt(sample.size-1)
    Gn.statistic[k] <- temp3 [k]^2 / Vk [k] 
  } 
  
  Choi.Shin.statistic <- max(Gn.statistic) 
  return ( Choi.Shin.statistic )
}


################################################################################
# mean of each column 
ColMean.fcn <-function(input.matrix=minus.thetaN.series, t1=1, t2= k )
{
  temp1<- input.matrix[t1:t2, ] 
  if ( length( temp1) == dim(y)[2] ) { temp2 <- input.matrix[t1:t2, ] 
  } else { temp2 <- apply ( input.matrix[ t1:t2, ] , 2, mean ) }
  return(temp2) 
}


################################################################################
# generate theta.1.k with k ranging with 1 to sample size n 
theta.1.k.mean.fcn <- function( input.matrix = dcc.score) { 
  sample.size <- dim(input.matrix )[1] 
  temp <- input.matrix [1 , ] 
  for (ind in (2: sample.size) ) { 
    input.matrix.sub <- input.matrix [1: ind, ] 
    temp<- rbind( temp, apply(input.matrix.sub, 2, mean) ) 
  } 
  return(temp)
}


################################################################################
# Univariate GARCH simulation
GARCHSim <- function(omega, alpha, beta, sims) {
  #sims (integer): Number of simulations performed. 
  ## initialize the vector of simulated returns and variances
  ret = numeric(sims)
  sigma2 = numeric(sims)
  
  ## initialize the variance at time t = 1 with the unconditional variance.
  sigma2[1] = omega/(1.0 - alpha - beta)
  
  ## sample the first observations
  ret[1] = rnorm(1, mean = 0, sd = sqrt(sigma2[1]))
  
  ## loop over sims.
  for (t in 2:sims) {
    #update volatility
    sigma2[t] = omega + alpha * ret[t - 1]^2 + beta * sigma2[t - 1]
    
    # sample new observarion
    ret[t] = rnorm(1, mean = 0, sd = sqrt(sigma2[t]))
  }
  
  ## we return a list with three components: the sampled returns, the variance, 
  ## and the standardized residuals. 
  lOut = list()
  lOut[["ret"]] = ret
  lOut[["sigma2"]] = sigma2
  
  # gives standardized residuals for demeaned returns. 
  lOut[["stdres"]] = ret/sqrt(sigma2)
  
  return(lOut)
} 


#y <- GARCHmultisim (params =rbind( t (as.matrix(c (0.1, 0.08, 0.89 ))), t (as.matrix(c (0.02, 0.05, 0.94 )))), sims=250 )
################################################################################
# Multiple simulations of the GARCH model.
GARCHmultisim <- function(params, sims){
  multisim_res <- matrix(0L, nrow=sims, ncol = nrow(params))
  multisim_sigma <- matrix(0L, nrow=sims, ncol = nrow(params))
  # rows in parameter matrix determines amount of garch simulations.
  for(i in 1:nrow(params)){
    # matrix() transforms it into a vector. 
    multisim_res[, i] <- matrix(GARCHSim(params[i,1], params[i,2], params[i,3], sims)[["stdres"]])
    multisim_sigma[, i] <- matrix(GARCHSim(params[i,1], params[i,2], params[i,3], sims)[["sigma2"]])
  }
  lOut = list()
  lOut[["res"]] = multisim_res
  lOut[["sigma"]] = multisim_sigma
  return(lOut)
}

 
################################################################################
#y3 <- DCCsim ( GARCH_params= rbind( t (as.matrix(c (0.1, 0.08, 0.89 ))), t (as.matrix(c (0.02, 0.05, 0.94 )))), a = 0.1, b=0.5, sims=250)
# Simulation for the Dynamic Conditional Correlations (DCC).
DCCsim <- function(GARCH_params, a, b, sims){
  stdres <- GARCHmultisim(GARCH_params, sims)[["res"]]
  Gsigs <- GARCHmultisim(GARCH_params, sims)[["sigma"]]
  
  mQ = cor(stdres)
  iN = ncol(stdres)
  iT = nrow(stdres)
  
  # initialize the array for the correlations
  aCor = array(0, dim = c(iN, iN, iT))
  
  # initialize the array for the Q matrices
  aQ = array(0, dim = c(iN, iN, iT))
  
  ## initialization at the unconditional cor
  aCor[,, 1] = mQ
  aQ[,,1] = mQ
  
  # main loop
  for (t in 2:iT) {
    # update the Q matrix
    aQ[,, t] = mQ * (1 - a - b) + a * t(stdres[t - 1, , drop = FALSE]) %*%
      stdres[t - 1, , drop = FALSE] + b * aQ[,,t - 1]
    
    ## Compute the correlation as Q_tilde^{-1/2} Q Q_tilde^{-1/2} 
    aCor[,, t] = diag(sqrt(1/diag(aQ[,, t]))) %*% aQ[,, t] %*% 
      diag(sqrt(1/diag(aQ[,, t]))) 
  }
  
  #Getting univariate variance matrix D_t:
  msigma <- array(0L, dim = c(ncol(Gsigs), ncol(Gsigs), nrow(Gsigs)))
  
  for(i in 1:ncol(Gsigs)){
    msigma[i, i, ] <- Gsigs[, i]
  }
  
  # constructing the covariance matrix
  DCC_COVAR <- array(0, dim = c(ncol(Gsigs), ncol(Gsigs), nrow(Gsigs)))
  DCC_returns <- matrix(0, nrow = nrow(Gsigs), ncol = ncol(Gsigs))
  
  for(i in 1:sims){
    DCC_COVAR[,,i] <- msigma[,,i]^0.5 %*% aCor[,,i] %*% msigma[,,i]^0.5
    # Getting returns:
    DCC_returns[i,] <- mvrnorm(1, rep(0,ncol(Gsigs)), DCC_COVAR[,,i])
  }
  
  lOut = list()
  
  # correlations
  lOut[["aCor"]] = aCor
  
  # univariate variances
  lOut[["Gsigs"]] = Gsigs
  
  # covariances
  lOut[["DCC_COVAR"]] = DCC_COVAR
  
  # returns
  lOut[["DCC_returns"]] = DCC_returns
  
  return(lOut)
}
