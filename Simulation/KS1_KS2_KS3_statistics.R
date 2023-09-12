
# Similar to Kiefer, Vogelsang and
# Bunzel (2000) and Lobato (2001), the Brownian motion is approximated by normalized
# sum of a sequence of i.i.d increments using 200000 N(0; 1) variables. 


# KS1 - the KS test statistic using small-b asymptotics and Bartlett kernel
# KS2 - the KS statistic based on naive self-normalization wt 
# KS3 - the adjusted-range based KS test statistic

# New revision to include various kernels and considered bandwidth selection both of Andrews and NW

################################################################################
rm(list=ls())  
set.seed(123456789 )
folder <- c("/Users/sunjiajing/R2023/Codes-simulation-study")
library(stats)


numrep<- 10000 

KS_range.loop <- rep(0, numrep)
KS_wt.loop <- rep(0, numrep)
KS_small_b.loop <- rep(0, numrep) 

for(i in 1:numrep)
{ 
  
  sample.size <- 200000
  times <- seq(0, 1, length.out=sample.size)
  
  sample.size <-length(times) 
  dW <- rnorm(sample.size)/sqrt(sample.size)
  W <- cumsum(dW)
  Brownian.Bridge <- W - times * W[sample.size]  #  The Brownian bridge from (0,0) to (1,target) 
  
  KS_small_b.loop[i]  <- max(abs(Brownian.Bridge)) 
  
  temp1 <- max(Brownian.Bridge)-min(Brownian.Bridge)
  KS1 <- max(abs(Brownian.Bridge))/ temp1 
  KS_range.loop[i] <-KS1 
  
  dt <- times[2]-times[1]
  
  Brownian.Bridge.integral <- 0 
  for (j in 1:length(Brownian.Bridge))
  {
    Brownian.Bridge.integral  <-   Brownian.Bridge.integral+ Brownian.Bridge[j]^2 *  dt ;
  }
  
  KS3 <- max(abs(Brownian.Bridge))/  sqrt(Brownian.Bridge.integral ) 
  KS_wt.loop[i]<-KS3 
  
  
  cat(i,"th replication" ,"\n") 
  
}


write.csv(KS_range.loop ,"KS_range_results.csv")   
write.csv(KS_small_b.loop, "KS_small_b_results.csv")  
write.csv(KS_wt.loop, "KS_wt_results.csv")  

#########################  generate histogram    ###########################   
jpeg('hist_small_b.png')
hist(KS_small_b.loop, 
     main= expression(), 
     xlab= expression(), 
     border="black", 
     col="grey", 
     breaks=20, prob=TRUE) # prob=TRUE for probabilities not counts

lines(density(KS_small_b.loop, adjust=2), col="blue", lwd=2)  
dev.off()

# main=expression('Histogram of Sup'[s]*'|B(s)|' ), 

#########################  generate histogram    ###########################   
jpeg('hist_KS_wt.png')
hist(KS_wt.loop, 
     main= expression(), 
     xlab= expression(), 
     col="grey", 
     breaks=20, prob=TRUE) 

lines(density(KS_wt.loop, adjust=2), col="blue", lwd=2)    

dev.off()

#########################  generate histogram    ###########################   

#library(EnvStats) 

# The adjusted-range based KS test statistic is bounded by 1, "density" uses the kernel based estimation method, which suffers
# from boundary value problem. Thus, here we fitted a smooth curve using the "loess" function. 

KS_range.loop<- read.csv("KS_range_results.csv") 
KS_range.loop <-  KS_range.loop[,2] 

temp <-hist(KS_range.loop, 
            main= expression(), 
            xlab= expression(), 
            border="black", 
            col="grey", 
            breaks=100, prob=TRUE) 
dev.off() 
y <- temp$density 
x <- 1: length(y) 
values <- loess(y ~ x) # fit a smooth curve for (x,y) 

jpeg('hist_KS_range-new.png')
hist(KS_range.loop, 
     main= expression(), 
     xlab= expression(), 
     border="black", 
     col="grey", 
     breaks=20, prob=TRUE) # prob=TRUE for probabilities not counts y <- temp$density  

breaks <- seq ( from = min(KS_range.loop ), to = max(KS_range.loop ), length.out = length(y))   

lines ( breaks, values$fitted , type="l", col="blue", lwd=2)    

dev.off() 