# date of coding: 5th of April 2023 
# suitability of the constant correlation (CC) model 

rm(list=ls())
set.seed(123456789 ) 
################################################################### 
folder <- c("/Users/sunjiajing/R2023/Codes-empirical-analysis")
setwd(folder)
source("source-empirical.R")

library(readxl)
library(forecast)
library(timeDate)
library(tseries)
library(sandwich)
library(xtable)

temp.dat <- read_xlsx ("data.xlsx") 
dat <- temp.dat [ setdiff ((1: dim(temp.dat)[1]),which (temp.dat $weekdays=="Saturday" | temp.dat $weekdays=="Sunday") ) , ]

Dow_Jones <- dat$Dow_Jones_Index_Industrial_Average
SP<- dat$`S&P/TSX_Composite_Index` 
FTSE_100 <- dat$Index_FTSE_100
CAC_40 <- dat$France_Index_CAC_40
DAX<- dat$Germany_Index_DAX 

Dow_Jones.return <- 	price2return(	dat$Dow_Jones_Index_Industrial_Average)
SP.return<- 	price2return(	dat$`S&P/TSX_Composite_Index` )
FTSE_100.return <- 	price2return(	dat$Index_FTSE_100)
CAC_40.return <- 	price2return(	dat$France_Index_CAC_40)
DAX.return<- 	price2return(	dat$Germany_Index_DAX )
################################################################################ 
temp <- cbind( Dow_Jones.return ,
               SP.return,
               FTSE_100.return ,
               CAC_40.return ,
               DAX.return )

loop.length <- 500
################################################################################
library(stats)
library(doParallel)
library(doSNOW)

no.cores<- detectCores() 
#numrep<- 2
numrep<- (dim(temp )[1]- loop.length+1) 
cl <- parallel::makeCluster(no.cores -1)
registerDoSNOW(cl) 
iterations <- numrep
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress) 

################################################################################ 
 result <- foreach (i = 1:(dim(temp )[1]- loop.length+1) , .combine = 'rbind' ,.options.snow = opts ) %dopar% 
{ 
    source("source-empirical.R")
 
Dow_Jones.return.interval   <- Dow_Jones.return[i: (i+loop.length-1)] 
SP.return.interval  <-SP.return [i: (i+loop.length-1)]
FTSE_100.return.interval  <-   FTSE_100.return[i: (i+loop.length-1)]
CAC_40.return.interval   <- CAC_40.return[i: (i+loop.length-1)]
DAX.return.interval   <- DAX.return [i: (i+loop.length-1)]

Dow_Jones.standardized.return.interval   <-  standardized.r  (y=    Dow_Jones.return.interval )  
SP.standardized.return.interval   <- standardized.r   (y=    SP.return.interval  ) 
FTSE_100.standardized.return.interval <-   standardized.r   (y=   FTSE_100.return.interval  ) 
CAC_40.standardized.return.interval   <-  standardized.r   (y= CAC_40.return.interval   ) 
DAX.standardized.return.interval   <-  standardized.r   (y=     DAX.return.interval  )   

################################################################################
# correlation coefficients. 

a.t<- cbind(    Dow_Jones.standardized.return.interval  ,  
                SP.standardized.return.interval, 
                FTSE_100.standardized.return.interval ,
                CAC_40.standardized.return.interval ,  
                DAX.standardized.return.interval       ) 

sample.size <- dim(  a.t) [1 ] ;  dim.no <- dim(  a.t) [2 ]  
temp.ind <- upper.tri ( diag(dim.no) )  

################################################################################ 
# G test of Shao and Zhang (2010) 
cor.vector <- c()

sample.size <- dim( a.t)[1]

for (cor.index in  (3: sample.size)) {  cor.vector <- rbind(cor.vector,     cor(a.t[1:cor.index,]  )  [   temp.ind])}




T.n.k.tilda   <- apply(  cor.vector  ,2, T.n.k.tilda.fcn )   # apply a function on each column of a matrix  

library(MASS)

tn.k <-  1/sqrt(dim(cor.vector )[1])*   T.n.k.tilda  


Gt.vec <- c( )

for ( k in 1: (dim( cor.vector)[1]-1))
{ 
  
  vn.k.component.1 <-  matrix( rep(0,  (dim( cor.vector)[2]* dim( cor.vector)[2])   ), nrow= dim( cor.vector)[2]) 
  
  for ( j in 1:k ) { 
    temp1 <-       corr.mat.fcn ( input.matrix =  a.t, t1=1, t2= j     ) -    corr.mat.fcn ( input.matrix =  a.t, t1=1, t2= k      ) 
    temp2 <- temp1 %*% t(temp1)  * j^2   
    vn.k.component.1 <- vn.k.component.1 +  temp2
  } 
  
  
  vn.k.component.2 <-  matrix( rep(0,  (dim( cor.vector)[2]* dim( cor.vector)[2])   ), nrow= dim( cor.vector)[2]) 
  
  for ( j in (k+1):   dim( cor.vector)[1]   ) {
    temp1 <-       corr.mat.fcn ( input.matrix =  a.t, t1= j  , t2=   dim( cor.vector)[1]     ) -    corr.mat.fcn ( input.matrix =  a.t, t1= (k+1), t2=   dim( cor.vector)[1]        ) 
    temp2 <- temp1 %*% t(temp1) *  (dim( cor.vector)[1]-j+1)^2 
    vn.k.component.2 <- vn.k.component.2 +  temp2
  } 
  
  v<- (vn.k.component.1+ vn.k.component.2)  /(dim( cor.vector)[1]^2)
  Gt  <-     tn.k   %*% ginv(v) %*%  t ( tn.k   ) 
  Gt.vec <- c(Gt.vec, Gt)
} 

Gt.correlation.matrix.statistic   <- max(  Gt.vec ) 

################################################################################ 
# correlation matrix test for choi and shin (2020)
Choi.Shin.correlation.matrix.statistic  <-  choi.shin.cor.fcn ( cor.vector =cor.vector )

################################################################################ 
# moderate the dependence using the ldl decomposition of sample variance of a.t 

sample.omega <- var(a.t) 
library(fastmatrix) 
temp <- ldl(sample.omega)
A.estimate<- temp$lower 
a.t <- t( ginv(A.estimate ) %*% t( a.t)) 

# generate correlation matrix 
cor.vector <- c()
for (cor.index in  (3: sample.size)) {  cor.vector <- rbind(cor.vector,     cor(a.t[1:cor.index,]  )  [   temp.ind])}

# adjusted-range based EKS test
T.n.k.tilda   <- apply(  cor.vector  ,2, T.n.k.tilda.fcn )   # apply a function on each column of a matrix 
v  <-    ( apply(  T.n.k.tilda, 2, max) - apply ( T.n.k.tilda , 2, min) )   
v.matrix <- diag(v^(-2))  
EKS.range.matrix  <- max (( T.n.k.tilda  ) %*% v.matrix  %*%  t (T.n.k.tilda  ) ) 

temp<- c( EKS.range.matrix , 
          Gt.correlation.matrix.statistic , 
          Choi.Shin.correlation.matrix.statistic )    
  } 

close(pb)
stopCluster(cl) 

################################################################################
folder.new <- paste0(folder, "/output")
setwd(folder.new )

temp.name <- paste0( "results_CC.csv") 
write.csv(result, temp.name) 



