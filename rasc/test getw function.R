get.w=function(nlk,psi,ncomm,nloc){
  soma=matrix(NA,nloc,ncomm-1)
  ge=ngreater(nlk=nlk,nloc=nloc,ncommun=ncomm) #greater or equal
  pnorm1=pnorm(0,mean=psi,sd=1)
  
  soma=getw(ge=ge, psi=psi, pnorm1=pnorm1,nlk=nlk,nloc=nloc,ncomm=ncomm)
  # for (i in 1:nloc){
  #   for (k in 1:(ncomm-1)){
  #     #for zil>k
  #     gt.k=ge[i,k+1]
  #     w.gt.k=tnorm.neg(n=gt.k,mu=psi[i,k],sig=1,pnorm1=pnorm1[i,k])
  #     
  #     #for zil=k
  #     eq.k=nlk[i,k]
  #     w.eq.k=tnorm.pos(n=eq.k,mu=psi[i,k],sig=1,pnorm1=pnorm1[i,k])
  #     
  #     #summarize
  #     soma[i,k]=sum(c(w.gt.k,w.eq.k))
  #   }
  # }
  list(soma=soma,nge=ge)
}    

get.w1=function(nlk,psi,ncomm,nloc){
  soma=matrix(NA,nloc,ncomm-1)
  ge=ngreater(nlk=nlk,nloc=nloc,ncommun=ncomm) #greater or equal
  pnorm1=pnorm(0,mean=psi,sd=1)
  
  for (i in 1:nloc){
    for (k in 1:(ncomm-1)){
      #for zil>k
      gt.k=ge[i,k+1]
      w.gt.k=tnorm.neg(n=gt.k,mu=psi[i,k],sig=1,pnorm1=pnorm1[i,k])

      #for zil=k
      eq.k=nlk[i,k]
      w.eq.k=tnorm.pos(n=eq.k,mu=psi[i,k],sig=1,pnorm1=pnorm1[i,k])

      #summarize
      soma[i,k]=sum(c(w.gt.k,w.eq.k))
    }
  }
  list(soma=soma,nge=ge)
}    
tnorm.pos <- function(n,mu,sig,pnorm1){   #generates truncated normal variates based on cumulative normal distribution
  z <- runif(n,pnorm1,1)
  qnorm(z,mu,sig)
}
tnorm.neg <- function(n,mu,sig,pnorm1){   #generates truncated normal variates based on cumulative normal distribution
  z <- runif(n,0,pnorm1)
  qnorm(z,mu,sig)
}

set.seed(1)
rcpp.res=get.w(nlk=nlk,psi=psi,ncomm=ncomm,nloc=nloc)
set.seed(1)
r.res=get.w1(nlk=nlk,psi=psi,ncomm=ncomm,nloc=nloc)
sum(round(rcpp.res$soma,4)!=round(r.res$soma,4))


