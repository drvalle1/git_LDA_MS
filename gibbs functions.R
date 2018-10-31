#this function generates truncated normal random variables
tnorm <- function(n,lo,hi,mu,sig){   
  #normal truncated lo and hi
  q1 <- pnorm(lo,mu,sig) #cumulative distribution
  q2 <- pnorm(hi,mu,sig) #cumulative distribution
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  z[z == -Inf]  <- lo
  z[z == Inf]   <- hi
  z
}
#-----------------------------------------
#this function generates vmat, which is then used to generate the theta matrix
get.delta.theta=function(nlk,gamma,ncomm,nloc,delta,sig2,mu){
  ngreater1=ngreater(nlk,nloc,ncomm)
  
  for (i in 1:nloc){
    for (j in 1:(ncomm-1)){
      #generate omega's
      omega.pos=tnorm(nlk[i,j]      ,lo=0,hi=10   ,mu=delta[i,j],sig=1)
      omega.neg=tnorm(ngreater1[i,j+1],lo=-10,hi=0,mu=delta[i,j],sig=1)
      
      #generate deltas
      prec=ngreater1[i,j+1]+(1/sig2)
      var1=1/prec
      soma=sum(c(omega.pos,omega.neg))
      pmedia=soma+(mu/sig2)
      delta[i,j]=rnorm(1,mean=var1*pmedia,sd=sqrt(var1))
    }
  }
  vmat=cbind(pnorm(delta),1)
  theta=convertVtoTheta(vmat,rep(1,nloc))
  list(theta=theta,delta=delta)  
}

#this function generates dirichlet random variables (1 one for each row of alpha)
rdirichlet1=function(alpha,ncomm,nspp){
  tmp=matrix(rgamma(n=ncomm*nspp,alpha,1),ncomm,nspp)
  soma=matrix(rowSums(tmp),ncomm,nspp)
  tmp/soma
}
