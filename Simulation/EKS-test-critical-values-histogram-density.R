# This code simulates the adjusted-range based KS and EKS test statistics and generates the critical values, 
# and plots the histograms and empirical density functions. 
################################################################################

rm(list=ls()) 
set.seed(123456789 )
folder <- c("/Users/sunjiajing/R2023/Codes-simulation-study")
setwd(folder)

numrep<- 10000
sample.size <- 5000

################################################################################
library(stats)
library(doParallel)
library(doSNOW) 
no.cores<- detectCores() 
cl <- makeCluster(no.cores -1)
registerDoSNOW(cl)
iterations <- numrep 
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress) 

result <- foreach (i = 1: numrep, .combine = 'rbind',.options.snow = opts ) %dopar% 
  { 
    sample.size <- length(times <- seq(0, 1, length.out = sample.size))
    
    generate_bridge <- function(){
      dW <- rnorm(sample.size) / sqrt(sample.size)
      W <- cumsum(dW)
      return(W - times * W[sample.size])
    }
    
    # Generate 10 Brownian Bridges
    num_bridges <- 10
    brownian.bridges <- replicate(num_bridges, generate_bridge())
    
    # Compute the self-normalizers
    v <- apply(brownian.bridges, 2, function(bridge) max(bridge) - min(bridge))
    
    # Convert to matrix
    Brownian.Bridge.matrix <- t(brownian.bridges) 
    
    v.matrix <- diag(v^(-2)) 
    
    ks.n.R <- numeric(num_bridges-1)
    
    for (i in 2:num_bridges  ) {
      ks.n.R[i-1] <- max(t(Brownian.Bridge.matrix[1:i, ]) %*% v.matrix[1:i, 1:i] %*% Brownian.Bridge.matrix[1:i, ])
    }
    
    temp <- ks.n.R
    
    
  }  

close(pb)
stopCluster(cl) 
################################################################################
folder <- paste0(folder, "/output")
setwd(folder)
temp_filename <- paste0( "ks_n_R_results_summary_q_1_10_numrep_", numrep, "_sample_size_", sample.size, ".csv" )
write.csv(result, temp_filename)


################################################################################
# histograms and empirical density function 

temp <- read.csv(temp_filename)
result <-temp[, 2:10]
row.names <- c( "90%", "95%", "97.5%", "99%", "99.5%", "99.9%") 
col.names <- paste0(" KS range q = ", 2:10) 

level.list <- c(0.900, 0.950, 0.975, 0.990, 0.995, 0.999)
temp.latex <- c()
for ( q in 1:9 ) 
{ 
  temp <- sort ( result [,q])
  temp.latex <- cbind(temp.latex, temp [level.list*numrep] )
}
rownames(temp.latex) <- row.names
colnames(temp.latex) <- col.names
temp.latex <- t(temp.latex) 

# generates latex output 
library(xtable)
temp_filename <- paste0( "ks_n_R_results_summary_q_1_10.tex" )
table1 <- xtable(temp.latex, digits= 4) 
print(table1 , file = temp_filename) 

################################################################################
# histograms and empirical density functions

for ( q in 2:10)
{ 
  tempfilename <- paste0( "ks_n_R_q_", q , "_numrep_", numrep, "_sample_size_", sample.size, ".png" )
  jpeg(tempfilename) 
  
  hist(result [,q] , 
       main= expression(), 
       xlab= expression( ), 
       ylab= expression( ), 
       border="black", 
       col="grey", 
       breaks=20, prob=TRUE) # prob=TRUE for probabilities not counts
  
  lines(density(result [,q] , adjust=2), col="blue", lwd=2) 
  dev.off()
}



