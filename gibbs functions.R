get.theta.from.lambda=function(lambda,ncomm){
  theta=rep(NA,ncomm)
  theta[1]=lambda[1]
  for (i in 2:ncomm){
    theta[i]=lambda[i]/prod(theta[1:(i-1)])
  }
  theta
}
#---------------------------------------------
#this function generates dirichlet random variables (1 one for each row of alpha)
rdirichlet1=function(alpha,ncomm,nspp){
  tmp=matrix(rgamma(n=ncomm*nspp,alpha,1),ncomm,nspp)
  soma=matrix(rowSums(tmp),ncomm,nspp)
  tmp/soma
}
#---------------------------------------------
ldirichlet=function(x,alpha){
  n=ncol(alpha)
  tmp=rowSums((alpha-1)*log(x))
  invBeta1=lgamma(alpha*n)-n*lgamma(alpha)
  invBeta1+tmp
}
#---------------------------------------------
ldmultinom1=function(size,x,prob){
  lgamma(size + 1) + sum(x * log(prob) - lgamma(x + 1))
}
#---------------------------------------------
get.media=function(theta,xmat,betas,ncomm,nloc){
  pmedia=exp(xmat%*%betas)
  lambda=GetLambda(LogTheta=log(theta),ncomm=ncomm)
  lambda1=matrix(lambda,nloc,ncomm,byrow=T)
  lambda1*pmedia
}
#---------------------------------------------
sample.theta=function(theta,a1,a2,b1,b2,nlk,ncomm,nloc,xmat,betas){
  for (i in 1:ncomm){
    soma.n=sum(nlk[,i:ncomm])
    medias=get.media(theta,xmat,betas,ncomm,nloc)
    p2=sum(medias[,i:ncomm])/theta[i]
    if (i==1) theta[i]=rgamma(1,soma.n+a1,p2+b1)
    if (i> 1) theta[i]=rgamma(1,soma.n+a2,p2+b2)
    if (theta[i]<0.0000001) theta[i]=0.0000001 #to avoid numerical problems
  }
  theta
}
#----------------------------------------------
sample.betas=function(lambda,nlk,xmat,betas,ncomm,nparam,jump1){
  betas.orig=betas.old=betas.prop=betas
  betas.prop[]=rnorm(nparam*ncomm,mean=betas.old,sd=jump1)
  for (i in 1:nparam){
    betas.new=betas.old
    betas.new[i,]=betas.prop[i,]
    
    lmedia.old=xmat%*%betas.old
    lmedia.new=xmat%*%betas.new
    p1.old=apply(nlk*lmedia.old,2,sum)
    p1.new=apply(nlk*lmedia.new,2,sum)
    
    p2.old=-lambda*apply(exp(lmedia.old),2,sum)
    p2.new=-lambda*apply(exp(lmedia.new),2,sum)

    p3.old=(1/2)*(betas.old[i,]^2)
    p3.new=(1/2)*(betas.new[i,]^2)
    
    pold=p1.old+p2.old+p3.old
    pnew=p1.new+p2.new+p3.new
    k=acceptMH(pold,pnew,betas.old[i,],betas.new[i,],F)
    betas.old[i,]=k$x
  }
  
  list(betas=betas.old,accept=betas.old!=betas.orig)
}

#---------------------------------
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
print.adapt = function(accept1z,jump1z,accept.output){
  accept1=accept1z; jump1=jump1z; 
  
  for (k in 1:length(accept1)){
    z=accept1[[k]]/accept.output
    print(names(accept1)[k])
    print(mean(z)); print(mean(jump1[[k]]))
  }
  
  for (k in 1:length(jump1)){
    cond=(accept1[[k]]/accept.output)>0.4 & jump1[[k]]<100
    jump1[[k]][cond] = jump1[[k]][cond]*2       
    cond=(accept1[[k]]/accept.output)<0.2 & jump1[[k]]>0.001
    jump1[[k]][cond] = jump1[[k]][cond]*0.5
    accept1[[k]][]=0
  }
  
  return(list(jump1=jump1,accept1=accept1))
}