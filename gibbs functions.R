fix.MH=function(lo,hi,old1,new1,jump){
  jold=pnorm(hi,mean=old1,sd=jump)-pnorm(lo,mean=old1,sd=jump)
  jnew=pnorm(hi,mean=new1,sd=jump)-pnorm(lo,mean=new1,sd=jump)
  log(jold)-log(jnew) #add this to pnew
}
#-------------------------------
#draw n samples from a multivariate normal distribution with mean 0 and covariance matrix sigma
#- n: number of samples to be generated
#- Sigma: covariance matrix
rmvnorm1=function (n, sigma, pre0.9_9994 = FALSE) 
{
  s. <- svd(sigma)
  if (!all(s.$d >= -sqrt(.Machine$double.eps) * abs(s.$d[1]))) {
    warning("sigma is numerically not positive definite")
  }
  R = t(s.$v %*% (t(s.$u) * sqrt(s.$d)))
  retval <- matrix(rnorm(n * ncol(sigma)), nrow = n, byrow = !pre0.9_9994) %*% R
  retval
}
#-------------------------------
tnorm <- function(n,lo,hi,mu,sig){   #generates truncated normal variates based on cumulative normal distribution
  #normal truncated lo and hi
  
  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))
  
  q1 <- pnorm(lo,mu,sig) #cumulative distribution
  q2 <- pnorm(hi,mu,sig) #cumulative distribution
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  z[z == -Inf]  <- lo[z == -Inf]
  z[z == Inf]   <- hi[z == Inf]
  z
}
#----------------------------------------------------------------------------------------------
acceptMH <- function(p0,p1,x0,x1,BLOCK){   #accept for M, M-H
  # if BLOCK, then accept as a block,
  # otherwise, accept individually
  
  nz           <- length(x0)  #no. to accept
  if(BLOCK) nz <- 1
  
  a    <- exp(p1 - p0)       #acceptance PR
  z    <- runif(nz,0,1)
  keep <- which(z < a)
  
  if(BLOCK & length(keep) > 0) x0 <- x1
  if(!BLOCK)                   x0[keep] <- x1[keep]           
  accept <- length(keep)        
  
  list(x = x0, accept = accept)
}
#----------------------------
print.adapt = function(accept1,jump1,accept.output){
  for (k in 1:length(accept1)){
    z=accept1[[k]]/accept.output
    print(names(accept1)[k])
    print(mean(z)); print(mean(jump1[[k]]))
  }
  
  for (k in 1:length(jump1)){
    cond=(accept1[[k]]/accept.output)>0.4 & jump1[[k]]<10000
    jump1[[k]][cond] = jump1[[k]][cond]*2       
    cond=(accept1[[k]]/accept.output)<0.2 & jump1[[k]]>0.001
    jump1[[k]][cond] = jump1[[k]][cond]*0.5
    accept1[[k]][]=0
  }
  
  return(list(jump1=jump1,accept1=accept1))
}
#-------------------------------
sample.omega=function(omega,ncomm,nloc,jump,xmat,betas,nlk,prod,sig2,o.lo,o.hi){
  omega.orig=omega.old=omega
  #if omega is too low or high, it makes:
  #a) random walk MH very sticky (hard to get out of very large absolute numbers)
  #b) beta estimates meaningless
  proposed=matrix(tnorm(ncomm*nloc,lo=o.lo,hi=o.hi,mu=omega,sig=jump),nloc,ncomm)  
  media=xmat%*%betas
  k1=-1/(2*sig2)
  
  #add correction for truncated normal proposal
  fixMH=fix.MH(lo=o.lo,hi=o.hi,old1=omega.old,new1=proposed,jump=jump)
    
  for (i in 1:(ncomm-1)){
    omega.new=omega.old
    omega.new[,i]=proposed[,i]

    #get theta
    e.omega.old=exp(omega.old)
    vmat.old=e.omega.old/(1+e.omega.old)
    theta.old=convertVtoTheta(vmat.old,prod)
    e.omega.new=exp(omega.new)
    vmat.new=e.omega.new/(1+e.omega.new)
    theta.new=convertVtoTheta(vmat.new,prod)
    
    #calculate log-probabilityes    
    pold=rowSums(nlk*log(theta.old)) #multinomial log-probs
    pnew=rowSums(nlk*log(theta.new)) #multinomial log-probs
    prior.old=k1*(omega.old[,i]-media[,i]) #prior log-probs
    prior.new=k1*(omega.new[,i]-media[,i]) #prior log-probs

    #accept/reject MH
    k=acceptMH(pold+prior.old,pnew+prior.new+fixMH[,i],omega.old[,i],omega.new[,i],F)
    omega.old[,i]=k$x
  }
  list(omega=omega.old,accept=omega.orig!=omega.old)
}
#----------------------------------------
sample.betas=function(ncomm,omega,sig2,xtx,nparam){
  prec=(1/sig2)*xtx+diag(1,nparam)
  var1=solve(prec)
  mvnorm=rmvnorm1(ncomm-1, sigma=var1)
    
  pmedia=(1/sig2)*t.xmat%*%omega[,-ncomm]
  media=var1%*%pmedia
  t(mvnorm)+media
}
#-----------------------------------------
sample.sig2=function(xmat,betas,omega,sig2.b,a1){
  media=xmat%*%betas
  err=(omega[,-ncomm]-media)
  b1=(t(err)%*%err/2)+sig2.b
  1/rgamma(1,a1,b1)
}