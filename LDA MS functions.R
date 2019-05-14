#-------------------------------
#this function generates dirichlet random variables (1 one for each row of alpha)
rdirichlet1=function(alpha,ncomm,nspp){
  tmp=matrix(rgamma(n=ncomm*nspp,alpha,1),ncomm,nspp)
  soma=matrix(rowSums(tmp),ncomm,nspp)
  tmp/soma
}
#-----------------------------------------
#this function samples the regression slope parameters
get.betas=function(tx,xmat,ncomm,npar,soma,n,nloc){
  betas=matrix(NA,npar,ncomm-1)
  for (i in 1:(ncomm-1)){
    qtq=diag(n[,i])
    prec=tx%*%qtq%*%xmat+diag(1,npar)
    var1=solve(prec)
    pmedia=tx%*%soma[,i]
    betas[,i]=rmvnorm(1,var1%*%pmedia,var1)
  }
  
  media=xmat%*%betas
  vmat=cbind(pnorm(media),1)
  theta=convertVtoTheta(vmat,rep(1,nloc))
  list(betas=betas,theta=theta)
}
#-----------------------------------------
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
#----------------------------------------------------------------------------------------------
get.w=function(xmat,nlk,betas,ncomm,nloc){
  media=xmat%*%betas
  soma=matrix(NA,nloc,ncomm-1)
  ge=ngreater(nlk=nlk,nloc=nloc,ncommun=ncomm) #greater or equal
  for (i in 1:nloc){
    for (k in 1:(ncomm-1)){
      #for zil>k
      gt.k=ge[i,k+1]
      w.gt.k=tnorm(gt.k,lo=-Inf,hi=0,mu=media[i,k],sig=1)
      
      #for zil=k
      eq.k=nlk[i,k]
      w.eq.k=tnorm(eq.k,lo=0,hi=Inf,mu=media[i,k],sig=1)
      
      #summarize
      soma[i,k]=sum(c(w.gt.k,w.eq.k))
    }
  }
  list(soma=soma,n=ge)
}    
#-------------------------------------
