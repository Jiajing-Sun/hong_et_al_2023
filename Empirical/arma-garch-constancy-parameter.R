# This code considers the structural change of conditional mean equations for the ARMA(1,1)-GARCH model. 

################################################################################
rm(list=ls())  
set.seed(123456789 )  
folder <- c("/Users/sunjiajing/R2023/Codes-empirical-analysis")
setwd(folder)
library(truncnorm)
library(MASS)
library(readxl) 
library(forecast) 
library(timeDate)
library(tseries) 
library(ggplot2)
library(sandwich)
library(msm) 
library(foreach)

temp.dat <- read_xlsx ("data.xlsx")   
dat <- temp.dat [ setdiff ((1: dim(temp.dat)[1]),which (temp.dat $weekdays=="Saturday" | temp.dat $weekdays=="Sunday") ) , ]

################################################################################
source("source-empirical.R")

Dow_Jones <- dat$Dow_Jones_Index_Industrial_Average 
SP<- dat$`S&P/TSX_Composite_Index` 
FTSE_100 <- dat$Index_FTSE_100
CAC_40 <- dat$France_Index_CAC_40
DAX  <- dat$Germany_Index_DAX 

Dow_Jones.return <- 	price2return(	dat$Dow_Jones_Index_Industrial_Average) 
SP.return<- 	price2return(	dat$`S&P/TSX_Composite_Index` )
FTSE_100.return <- 	price2return(	dat$Index_FTSE_100)
CAC_40.return <- 	price2return(	dat$France_Index_CAC_40)
DAX.return  <- 	price2return(	dat$Germany_Index_DAX ) 

################################################################################

temp <- cbind( Dow_Jones.return ,
               SP.return,
               FTSE_100.return ,
               CAC_40.return ,
               DAX.return    )  

loop.length <- 500 

################################################################################  

library(stats)
library(doParallel)
library(doSNOW)

no.cores<- detectCores() 
numrep<- dim(temp )[2]* (dim(temp )[1]- loop.length+1)

cl <- makeCluster(no.cores -1)
registerDoSNOW(cl)
iterations <- numrep
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress) 


result  <- foreach (stock.ind   = 1: dim(temp )[2]   , .combine = 'cbind'  )  %:% 
  foreach (i = 1: (dim(temp )[1]- loop.length+1), .combine = 'rbind' ,.options.snow = opts  )  %dopar% 
  {    
    source("source-empirical.R")
    library(MASS) 
    library(tseries)
    library(fGarch)  
    y <- temp[,stock.ind] [i: (i+loop.length-1)]   
    arma.1.1<- arma(y, order=c(1,1))
    res<- na.omit (arma.1.1$res)
    
   # try  (fit <- garchFit(~ arma(1,1)+garch(1,1), data = y,include.mean=TRUE ) )
   # res<- try (fit@residuals )   
    
    sample.size <- length(res)
    average.res <-mean (res )  
    demeaned.res <-res - average.res 
    
    ################################################################################
   # adjusted range based KS 
    
    cumsum.range <- cumsum(demeaned.res)
    R.hat <- ( max(cumsum.range) - min (cumsum.range))/sqrt( sample.size)   
    KS.R <- (max(abs(cumsum.range))/sqrt( sample.size) )/R.hat  
    KS.R.acceptance.rejection.asymptotics  <- (  KS.R  <=  0.9117 ) 
    
    ################################################################################
    # KS based on the SN of Shao (2010) and Lobato (2001)
    
    KS.Ws.numerator <-   max(abs(cumsum.range))  /sqrt(sample.size) 
    temp1 <-  ( cumsum(res) / (1: length(res) )  -  average.res *rep(1, length(res)) ) ^2
    temp2 <- (1: length(res)) ^2
    temp3 <- temp1*temp2 
    wt.temp <- sum( temp3) /(sample.size^2)
    KS.W <-   KS.Ws.numerator  / sqrt( wt.temp )  
    KS.W.acception.rejection.asymptotics  <-(KS.W<=  3.0585   ) 
    
    ################################################################################
    # Gn    
    v.k <- rep( 0, sample.size-1) 
    Gn <- rep ( 0, sample.size-1) 
    for ( k in 1: (sample.size-1))
    { 
      temp1 <-0       
      for ( j in 1:k) { temp1 <- temp1+   (s_t1_t2.fcn(t1=1, t2=j  ) -j/k*  s_t1_t2.fcn(t1=1, t2=k ) ) ^2  }
      temp2 <-0 
      for ( j in (k+1):   (sample.size-1)   )  {   temp2  <- temp2 +  (s_t1_t2.fcn(t1=j, t2=  sample.size-1  ) - (sample.size - j+1)/(sample.size - k)*  s_t1_t2.fcn(t1= (k+1), t2=sample.size-1 ) ) ^2 } 
      v.k [k] <-    (temp1+ temp2)/    sample.size^2 
      demeaned.res.v2 <- demeaned.res [ 1: (sample.size-1)]
      temp3 <-  cumsum(demeaned.res.v2) /sqrt(sample.size-1)
      Gn[k] <-    temp3 [k]^2 /  v.k [k] 
    } 
    
    Gn  <- max(Gn) 
    Gn.acception.rejection.asymptotics <-     (Gn <=   40.1)
    
    ################################################################################
    # KS test statistics based on small-b (standard) asymptotics and Bartlett kernel 
    
    KS.0.statistics.numerator <-(max(abs(cumsum.range)) /sqrt(sample.size)) 
    library(cointReg) 
    s.e.Bartlett.denominator  <- getLongRunVar(res, kernel = "ba", bandwidth = "and")
    s.e.Bartlett.denominator  <- sqrt(s.e.Bartlett.denominator $Omega) 
    KS.0 <- KS.0.statistics.numerator/ s.e.Bartlett.denominator  
    KS.0.Bartlett.acception.rejection.asymptotics <-   (KS.0<= 1.3640  ) 
    
    temp1 <- c (KS.R,  
                Gn,
                KS.W, 
                KS.0 )  
    
    temp2 <- c (  KS.R.acceptance.rejection.asymptotics ,  
                  Gn.acception.rejection.asymptotics ,  
                  KS.W.acception.rejection.asymptotics, 
                  KS.0.Bartlett.acception.rejection.asymptotics  ) 
    
    temp3 <- c(temp1, temp2  )
  } 

close(pb)
stopCluster(cl) 
################################################################################
folder.new <- paste0(folder, "/output")
setwd(folder.new )

temp.name <- paste0( "results_parameter_constancy_conditional_mean_model_ARMA(1,1)-GARCH(1,1).csv") 
write.csv(result, temp.name) 


