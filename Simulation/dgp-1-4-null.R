# This codes generates the null for EKS and G under the following DGPs: 
# DGP1: Simple homoskedastic errors (two independent series stacked up). 
# DGP2: VAR with homoskedastic errors
# DGP3: VAR with conditional heteroskedastic errors
# DGP4: VAR with unconditional heteroskedastic errors

rm(list=ls()) 
set.seed(123456789 )
folder <- c("/Users/sunjiajing/R2023/Codes-simulation-study")
setwd(folder) 
start.time=proc.time()
################################################################################
numrep <- 1000
n <- 250 
dgp.type <- 4 # choose from 1 to 4 
################################################################################
start.time=proc.time()

library(stats)
library(doParallel)
no.cores= detectCores() 
cl =makeCluster(no.cores-1)
library(doSNOW)
registerDoSNOW(cl)

iterations = numrep 
pb = txtProgressBar(max = iterations, style = 3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress) 

################################################################################
result <- foreach (i = 1: numrep , .combine = 'rbind' ,.options.snow = opts ) %dopar% 
  { 
    ######################################### 
    # p: dimension of covariance matrix with parameter rho
    
    S<-function(rho,p)
    {
      a<-rho^(1:(p-1));
      b<-diag(1,p);
      for(i in 1:(p-1))
      {
        b[i,(i+1):p]<-a[1:(p-i)];
      }
      b<-b+t(b)-diag(1,p);
      b;
    } 
    ######################################### 
    sample.size <- n   
    library(MASS)
    
    if (dgp.type==1){ 
      Sigma <- diag(rep(1,2)) 
      alpha <- matrix(c(0.5, 0,0, 0.5), 2)
      e <- mvrnorm(n= sample.size,mu=rep(0,dimension.no),Sigma =Sigma) 
      y <- e
      for(j in 2: sample.size )
      {y[j,]<- alpha %*% y[j-1,]+ e[j,]} } 
    
    if (dgp.type==2){  
      rho <-0.1; dimension.no <- 2 
      Sigma <- S(rho,dimension.no) 
      alpha <- matrix(c(0.5, 0.1, 0.1, 0.5), 2) 
      e <- mvrnorm(n= sample.size,mu=rep(0,dimension.no),Sigma =Sigma) 
      y <- e
      for(j in 2: sample.size ) 
      {y[j,]<- alpha %*% y[j-1,]+ e[j,]} } 
    
    if (dgp.type==3){ 
      rho <-0.1; dimension.no <- 2 
      Sigma <- S(rho,dimension.no) 
      alpha <- matrix(c(0.5, 0.1, 0.1, 0.5), 2) 
      alpha1 <- 0.1 
      beta1 <- 0.89 
      e <- mvrnorm(n= sample.size,mu=rep(0,dimension.no),Sigma =Sigma) 
      sigma.t.squared <- 0*e
      sigma.t.squared [1, ] <- ( e [1, ])^2 
      epsilon.t <- 0*e
      for ( i in 2:( sample.size)) { 
        epsilon.t [i-1, ] <- sqrt ( sigma.t.squared[i-1, ] ) * e [i-1, ]
        sigma.t.squared [i, ] = (1 - alpha1-beta1) + alpha1 * ( epsilon.t [i-1, ] )^2 + beta1* sigma.t.squared [i-1, ] 
      } 
      e <- sqrt( sigma.t.squared) * e 
      y <- e
      for(j in 2: sample.size )
      {y[j,]<- alpha %*% y[j-1,]+ e[j,]} } 
    
    if (dgp.type==4) {
      rho <-0.1; dimension.no <- 2 
      Sigma <- S(rho,dimension.no) 
      alpha <- matrix(c(0.5, 0.1, 0.1, 0.5), 2) 
      delta <- 1 
      e <- mvrnorm(n= sample.size,mu=rep(0,dimension.no),Sigma =Sigma) 
      sigma.t.squared <- 0*e
      sigma.t.squared [1:sample.size/2,] <-1 
      sigma.t.squared [(1+sample.size/2):sample.size,] <-1 +delta 
      e <- sqrt( sigma.t.squared) * e 
      y <- e
      for(j in 2: sample.size )
      {y[j,]<- alpha %*% y[j-1,]+ e[j,]} 
    }
    
    demeaned.y <- y - colMeans(y)  
    
    ######################################### 
    # functions:
    # adjusted-range based EKS statistic
    ks.range.fcn <- function(input.vec = u.hat.cls.matrix) {
      colmeans.theta.1.k <- colMeans(input.vec ) 
      demeaned.input.vec <- input.vec - cbind( rep( colmeans.theta.1.k [1], sample.size) , 
                                               rep( colmeans.theta.1.k [2], sample.size) ) 
      cumsum.range <-apply( demeaned.input.vec , 2, cumsum)
      v<- diag (( apply( cumsum.range , 2, max) - apply ( cumsum.range , 2, min) )^(-2) )
      KS_range <- cumsum.range %*% v %*% t (cumsum.range )
      KS_range <- max( KS_range )}
    # mean of each column 
    ColMean.fcn <-function(input.matrix=minus.thetaN.series, t1=1, t2= k )
    {
      temp1<- input.matrix[t1:t2, ] 
      if ( length( temp1) == dim(y)[2] ) { temp2 <- input.matrix[t1:t2, ] 
      } else { temp2 <- apply ( input.matrix[ t1:t2, ] , 2, mean ) }
      return(temp2) 
    }
    # generate theta.1.k with k ranging with 1 to sample size n 
    theta.1.k.mean.fcn <- function( input.matrix = original.series) { 
      sample.size <- dim(input.matrix )[1] 
      temp <- input.matrix [1 , ] 
      for (ind in (2: sample.size) ) { 
        input.matrix.sub <- input.matrix [1: ind, ] 
        temp<- rbind( temp, apply(input.matrix.sub, 2, mean) ) 
      } 
      return(temp)
    }
    
    ######################################## 
    # G test statistic from Shao and Zhang (2010) 
    
    original.series <- y 
    theta.1.k <- theta.1.k.mean.fcn(input.matrix = original.series ) 
    tn.k<-(1:dim(original.series)[1])/sqrt(dim(original.series)[1])*(theta.1.k-cbind(rep(theta.1.k[dim(original.series)[1],1],dim(original.series)[1]),
                                                                                     rep(theta.1.k[dim(original.series)[1],2],dim(original.series)[1])))
    Gt.vec <- c( )
    for ( k in 1: (sample.size-1))
    { 
      vn.k.component.1 <- matrix( rep(0, (dim(original.series)[2]* dim(original.series)[2]) ), nrow= dim(original.series)[2]) 
      for ( j in 1:k ) { 
        temp1 <- ColMean.fcn ( input.matrix = original.series, t1=1, t2= j ) - ColMean.fcn ( input.matrix = original.series, t1=1, t2= k ) 
        temp2 <- temp1 %*% t(temp1) * j^2 
        vn.k.component.1 <- vn.k.component.1 + temp2
      } 
      vn.k.component.2 <- matrix( rep(0, (dim(original.series)[2]* dim(original.series)[2]) ), nrow= dim(original.series)[2]) 
      for ( j in (k+1): sample.size ) {
        temp1 <- ColMean.fcn ( input.matrix = original.series, t1= j , t2= sample.size ) - ColMean.fcn ( input.matrix = original.series, t1= (k+1), t2= sample.size ) 
        temp2 <- temp1 %*% t(temp1) * (sample.size-j+1)^2 
        vn.k.component.2 <- vn.k.component.2 + temp2
      } 
      v<- (vn.k.component.1+ vn.k.component.2) /(sample.size^2)
      Gt <- tn.k %*% ginv(v) %*% t ( tn.k ) 
      Gt.vec <- c(Gt.vec, Gt)
    } 
    Gt <- max( Gt.vec ) 
    
    ######################################### 
    # adjusted-range based EKS test statistics  
    #  use LDL decomposition of the sample variance 
    sample.omega <- var(y)
    library(fastmatrix) 
    temp <- ldl(sample.omega)
    A.estimate<- temp$lower 
    u.hat.cls.matrix <- t( ginv(A.estimate ) %*% t( y)) 
    sample.size <- dim( u.hat.cls.matrix )[1] 
    KS_range <- ks.range.fcn(input.vec=u.hat.cls.matrix)
    
    
    temp <- c( KS_range, Gt ) 
  } 

close(pb)
stopCluster(cl) 
################################################################################
# saving statistics values 

folder.new <- paste0( folder, "/dgp", dgp.type ) 
setwd(folder.new)
temp.name <- paste0( "dgp", dgp.type, "_null_sample_size_",n, "_num_rep_", numrep ,".csv") 
write.csv(result, temp.name) 
proc.time()-start.time 

