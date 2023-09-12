# This code tests the parameter stability of the bivariate DCC-GARCH(1,1) model

################################################################################

rm(list=ls())
set.seed(123456789 ) 
start_time <- Sys.time()

library(readxl)
library(forecast)
library(timeDate)
library(tseries)
library(sandwich)
library(xtable) 
library(imputeTS)
################################################################### 
folder <- c("/Users/sunjiajing/R2023/Codes-empirical-analysis")
setwd(folder) 
source("source-empirical.R") 

temp.dat <- read_xlsx ("data.xlsx") 
dat <- temp.dat [ setdiff ((1: dim(temp.dat)[1]),which (temp.dat $weekdays=="Saturday" | temp.dat $weekdays=="Sunday") ) , ] 

Dow_Jones <- dat$Dow_Jones_Index_Industrial_Average
NYSE <- dat$NYSE_Index_Composite
SP<- dat$S.P.TSX_Composite_Index
FTSE_100 <- dat$Index_FTSE_100
CAC_40 <- dat$France_Index_CAC_40
DAX<- dat$Germany_Index_DAX 
FTSE_MIB <- dat$Index_FTSE_MIB 
SSE <- dat$Shanghai_Stock_Exchange_Index_Composite
Nikkei_225 <- dat$Tokyo_Index_Nikkei_225
KOSPI_200 <- dat$KRX_Index_KOSPI_200 

Dow_Jones.return <- 	price2return(	dat$Dow_Jones_Index_Industrial_Average)
NYSE.return <- 	price2return(	dat$NYSE_Index_Composite)
SP.return<- 	price2return(	dat$S.P.TSX_Composite_Index )
FTSE_100.return <- 	price2return(	dat$Index_FTSE_100)
CAC_40.return <- 	price2return(	dat$France_Index_CAC_40)
DAX.return<- 	price2return(	dat$Germany_Index_DAX )
FTSE_MIB.return <- 	price2return(	dat$Index_FTSE_MIB)
SSE.return <- 	price2return(	dat$Shanghai_Stock_Exchange_Index_Composite)
Nikkei_225.return <- 	price2return(	dat$Tokyo_Index_Nikkei_225)
KOSPI_200.return <- 	price2return(	dat$KRX_Index_KOSPI_200) 

################################################################################
loop.length <- 250
numstep <- length(Dow_Jones.return) - loop.length+1 
# set.index <- seq(1,numstep,5) 
set.index <- 1:numstep 

file.name.string<-c("Dow_Jones.return", "SP.return", "FTSE_100.return", "CAC_40.return" , "DAX.return") 
combination.string <- combn(file.name.string, m=2)

library(cointReg) # estimate long run variance. 
library(stats)
library(doParallel)
library(doSNOW)

no.cores<- detectCores() 
cl <- makeCluster(no.cores -1)
registerDoSNOW(cl)
iterations <- length(set.index)* dim(combination.string )[2] 
# iterations <- length(numstep)* dim(combination.string )[2] 
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress) 


result <- foreach (string.ind = 1: 10 , .combine = 'cbind' ) %:%
  foreach (ind = 1: length (set.index), .combine = 'rbind' ,.options.snow = opts ) %dopar% 
  { 
    source( paste0(folder, "/source-empirical.R")) 
    i <- set.index [ind] 
    
    library(tseries)
    library(MASS) 
    
    if (string.ind==1) {interval.1 <- Dow_Jones.return [i: (i+loop.length-1)];interval.2 <- SP.return [i: (i+loop.length-1)] } 
    if (string.ind==2) {interval.1 <- Dow_Jones.return [i: (i+loop.length-1)];interval.2 <- FTSE_100.return [i: (i+loop.length-1)] } 
    if (string.ind==3) {interval.1 <- Dow_Jones.return [i: (i+loop.length-1)];interval.2 <- CAC_40.return [i: (i+loop.length-1)] } 
    if (string.ind==4) {interval.1 <- Dow_Jones.return [i: (i+loop.length-1)];interval.2 <- DAX.return [i: (i+loop.length-1)] } 
    if (string.ind==5) {interval.1 <- SP.return [i: (i+loop.length-1)];interval.2 <- FTSE_100.return[i: (i+loop.length-1)] } 
    if (string.ind==6) {interval.1 <- SP.return [i: (i+loop.length-1)];interval.2 <- CAC_40.return[i: (i+loop.length-1)] } 
    if (string.ind==7) {interval.1 <- SP.return [i: (i+loop.length-1)];interval.2 <- DAX.return[i: (i+loop.length-1)] } 
    if (string.ind==8) {interval.1 <- FTSE_100.return [i: (i+loop.length-1)];interval.2 <- CAC_40.return[i: (i+loop.length-1)] } 
    if (string.ind==9) {interval.1 <- FTSE_100.return [i: (i+loop.length-1)];interval.2 <- DAX.return[i: (i+loop.length-1)] } 
    if (string.ind==10) {interval.1 <- CAC_40.return [i: (i+loop.length-1)];interval.2 <- DAX.return[i: (i+loop.length-1)] } 
    
    results.1 <- arma.11.residuals (y= interval.1 ) 
    results.2 <- arma.11.residuals (y= interval.2 ) 
    returns.interval <- data.frame( results.1 , results.2 )
    
    ################################################################################
    # DCC estimation
    
    library(rugarch); library(rmgarch)
    
    garch11.spec = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = FALSE), 
                              variance.model = list(garchOrder = c(1,1), 
                                                    model = "sGARCH"), 
                              distribution.model = "norm") 
    
    # dcc specification - GARCH(1,1) for conditional correlations
    dcc.garch11.spec = dccspec(uspec = multispec( replicate(2, garch11.spec) ), 
                               dccOrder = c(1,1), 
                               distribution = "mvnorm")

    dcc.fit = dccfit(dcc.garch11.spec, data = returns.interval ) 
    # class(dcc.fit)
    # summary(dcc.fit)
    # slotNames(dcc.fit)
    # names(dcc.fit@mfit)
    # names(dcc.fit@model)
    # length(dcc.fit@mfit$coef)
    dcc.score <- dcc.fit@mfit$scores 
    colnames(dcc.score) <- paste0("y" , 1: dim(dcc.score)[2]) 
    y.store <- dcc.score
    
    
    ################################################################################
    # adjusted-range based constancy of parameter test 
    
    sample.omega <- var(dcc.score )
    temp <- fastmatrix::ldl(sample.omega) # use LDL decomposition of the sample variance of the prewhitened errors
    A.estimate<- temp$lower 
    u.hat.primitive.shocks.matrix <- t( ginv(A.estimate ) %*% t( dcc.score )) 
    sample.size <- dim( u.hat.primitive.shocks.matrix )[1] 
    
    theta.1.k <- theta.1.k.mean.fcn(input.matrix = u.hat.primitive.shocks.matrix ) 
    T.n.k.mean.tilda <- apply( theta.1.k ,2, T.n.k.tilda.fcn ) 
    v <- ( apply( T.n.k.mean.tilda , 2, max) - apply ( T.n.k.mean.tilda , 2, min) ) 
    v.matrix <- diag(v^(-2)) 
    EKS.range.matrix <- max (( T.n.k.mean.tilda ) %*% v.matrix %*% t (T.n.k.mean.tilda ) )
    
    ################################################################################
    # G test statistic of Shao and Zhang (2010) 
    
    y <- y.store 
    theta.1.k <- theta.1.k.mean.fcn(input.matrix = y ) 
    T.n.k.mean.tilda <- apply( theta.1.k ,2, T.n.k.tilda.fcn ) 
    tn.k<- T.n.k.mean.tilda / sqrt(sample.size) 
    Gt.vec <- c( ) 
    sample.size <- dim(y)[1]
    
    for ( k in 1: (sample.size-1))
    { 
      vn.k.component.1 <- matrix( rep(0, (dim(y)[2]* dim(y)[2]) ), nrow= dim(y)[2]) 
      for ( j in 1:k ) { 
        temp1 <- ColMean.fcn ( input.matrix = y, t1=1, t2= j ) - ColMean.fcn ( input.matrix = y, t1=1, t2= k ) 
        temp2 <- temp1 %*% t(temp1) * j^2 
        vn.k.component.1 <- vn.k.component.1 + temp2
      } 
      vn.k.component.2 <- matrix( rep(0, (dim(y)[2]* dim(y)[2]) ), nrow= dim(y)[2]) 
      for ( j in (k+1): sample.size ) {
        temp1 <- ColMean.fcn ( input.matrix = y, t1= j , t2= sample.size ) - ColMean.fcn ( input.matrix = y, t1= (k+1), t2= sample.size ) 
        temp2 <- temp1 %*% t(temp1) * (sample.size-j+1)^2 
        vn.k.component.2 <- vn.k.component.2 + temp2
      } 
      v<- (vn.k.component.1+ vn.k.component.2) /(sample.size^2)
      Gt <- tn.k %*% ginv(v) %*% t ( tn.k ) 
      Gt.vec <- c(Gt.vec, Gt)
    } 
    
    Gt <- max( Gt.vec ) 
    temp <- c(EKS.range.matrix, Gt) 
    
  } 


close(pb)
stopCluster(cl) 

################################################################################
folder.name <- paste0(folder , "/output"); setwd(folder.name) # this is used to store output!

end_time <- Sys.time()
running.time <- end_time - start_time
running.time 

#write.table(running.time, file = paste0("running-time_20221229.txt"))

temp.name <- paste0( "results_dcc-garch-", loop.length , ".csv") 
write.csv(result, temp.name) 
 