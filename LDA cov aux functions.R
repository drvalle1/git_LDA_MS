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
sample.betas=function(nlk,xmat,betas,ncomm,nparam,jump,a.gamma,b.gamma,var.betas){
  betas.orig=betas.old=betas.prop=betas
  betas.prop[]=rnorm(nparam*ncomm,mean=betas.old,sd=jump)
  p2=a.gamma+colSums(nlk)
  for (i in 1:nparam){
    betas.new=betas.old
    betas.new[i,]=betas.prop[i,]
    
    lmedia.old=xmat%*%betas.old
    lmedia.new=xmat%*%betas.new
    p1.old=colSums(nlk*lmedia.old)
    p1.new=colSums(nlk*lmedia.new)
    
    media.old=exp(lmedia.old)
    media.new=exp(lmedia.new)
    p3.old=log(b.gamma+colSums(media.old))
    p3.new=log(b.gamma+colSums(media.new))
    
    prior.old=(1/(2*var.betas))*(betas.old[i,]^2)
    prior.new=(1/(2*var.betas))*(betas.new[i,]^2)
    
    pold=p1.old-p2*p3.old+prior.old
    pnew=p1.new-p2*p3.new+prior.new
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
#----------------------------
get.lambda=function(nlk,a.gamma,b.gamma,xmat,betas,ncomm){
  nk=colSums(nlk)
  a1=nk+a.gamma
  media=exp(xmat%*%betas)
  b1=colSums(media)+b.gamma
  rgamma(ncomm,a1,b1)
}