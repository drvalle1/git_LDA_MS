get.w=function(nlk,psi,ncomm,nloc){
  soma=matrix(NA,nloc,ncomm-1)
  ge=ngreater(nlk=nlk,nloc=nloc,ncommun=ncomm) #greater or equal
  for (i in 1:nloc){
    for (k in 1:(ncomm-1)){
      #for zil>k
      gt.k=ge[i,k+1]
      w.gt.k=tnorm(n=gt.k,lo=-Inf,hi=0,mu=psi[i,k],sig=1)
      
      #for zil=k
      eq.k=nlk[i,k]
      w.eq.k=tnorm(n=eq.k,lo=0,hi=Inf,mu=psi[i,k],sig=1)
      
      #summarize
      soma[i,k]=sum(c(w.gt.k,w.eq.k))
    }
  }
  list(soma=soma,nge=ge)
}    
#-------------------------------------
tnorm <- function(n,lo,hi,mu,sig){   #generates truncated normal variates based on cumulative normal distribution
  #normal truncated lo and hi
  
  q1 <- pnorm(lo,mu,sig) #cumulative distribution
  q2 <- pnorm(hi,mu,sig) #cumulative distribution
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  z[z == -Inf]  <- lo
  z[z == Inf]   <- hi
  z
}
#----------------------------------
library(tictoc)

tic()
for (i in 1:10){
  get.w(nlk=nlk,psi=psi,ncomm=ncomm,nloc=nloc)  
}
toc()
#----------------------------


get.w1=function(nlk,psi,ncomm,nloc){
  soma=matrix(NA,nloc,ncomm-1)
  ge=ngreater(nlk=nlk,nloc=nloc,ncommun=ncomm) #greater or equal
  for (i in 1:nloc){
    for (k in 1:(ncomm-1)){
      #for zil>k
      gt.k=ge[i,k+1]
      w.gt.k=tnorm.neg(n=gt.k,mu=psi[i,k],sig=1)
      
      #for zil=k
      eq.k=nlk[i,k]
      w.eq.k=tnorm.pos(n=eq.k,mu=psi[i,k],sig=1)
      
      #summarize
      soma[i,k]=sum(c(w.gt.k,w.eq.k))
    }
  }
  list(soma=soma,nge=ge)
}    
#-------------------------------------
tnorm.pos <- function(n,mu,sig){   #generates truncated normal variates based on cumulative normal distribution
  #normal truncated lo and hi
  
  q1 <- pnorm(0,mu,sig) #cumulative distribution
  q2 <- 1
  
  z <- runif(n,q1,q2)
  qnorm(z,mu,sig)
}
#----------------------------------
tnorm.neg <- function(n,mu,sig){   #generates truncated normal variates based on cumulative normal distribution
  #normal truncated lo and hi
  
  q1 <- 0
  q2 <- pnorm(0,mu,sig) #cumulative distribution
  
  z <- runif(n,q1,q2)
  qnorm(z,mu,sig)
}
#----------------------------------

library(tictoc)

tic()
for (i in 1:10){
  get.w1(nlk=nlk,psi=psi,ncomm=ncomm,nloc=nloc)  
}
toc()
