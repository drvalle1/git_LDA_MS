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
#-------------------------------------
# tnorm.pos <- function(n,mu,sig,pnorm1){   #generates truncated normal variates based on cumulative normal distribution
#   z <- runif(n,pnorm1,1)
#   qnorm(z,mu,sig)
# }
# #----------------------------------
# tnorm.neg <- function(n,mu,sig,pnorm1){   #generates truncated normal variates based on cumulative normal distribution
#   z <- runif(n,0,pnorm1)
#   qnorm(z,mu,sig)
# }
#-------------------------------------
#this function generates dirichlet random variables (1 one for each row of alpha)
rdirichlet1=function(alpha,ncomm,nspp){
  tmp=matrix(rgamma(n=ncomm*nspp,alpha,1),ncomm,nspp)
  soma=matrix(rowSums(tmp),ncomm,nspp)
  tmp/soma
}
#---------------------------------------------
get.psi=function(soma,nge,sd0,mu0,ncomm){
  var0=sd0^2
  prec=nge[,-ncomm]+(1/var0)
  var1=1/prec
  pmedia=soma+(mu0/var0)
  tmp=rnorm(nloc*(ncomm-1),mean=var1*pmedia,sd=sqrt(var1))
  matrix(tmp,nloc,ncomm-1)
}