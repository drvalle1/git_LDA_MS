#---------------------------------------------
#this function generates dirichlet random variables (1 one for each row of alpha)
rdirichlet1=function(alpha,ncomm,nspp){
  tmp=matrix(rgamma(n=ncomm*nspp,alpha,1),ncomm,nspp)
  soma=matrix(rowSums(tmp),ncomm,nspp)
  tmp/soma
}
#---------------------------------------------
ldirichlet=function(x,alpha){
  n=ncol(x)
  tmp=rowSums((alpha-1)*log(x))
  invBeta1=lgamma(alpha*n)-n*lgamma(alpha)
  invBeta1+tmp
}
#---------------------------------------------
ldmultinom1=function(size,x,prob){
  lgamma(size + 1) + sum(x * log(prob) - lgamma(x + 1))
}
#---------------------------------------------
sample.betas=function(y,xmat,betas,ncomm,nparam,jump,var.betas,phi,ntot){
  betas.orig=betas.old=betas.prop=betas
  betas.prop[]=rnorm(nparam*ncomm,mean=betas.old,sd=jump)
  
  var.betas1=matrix(var.betas,nparam,ncomm)
  prior.old=dnorm(betas.orig,mean=0,sd=sqrt(var.betas1),log=T)
  prior.new=dnorm(betas.prop,mean=0,sd=sqrt(var.betas1),log=T)
  
  runif1=matrix(runif(nparam*ncomm),nparam,ncomm)
  betas=sampleBetas(y=y,xmat=xmat,betas_prop=betas.prop,
                    prior_old=prior.old, prior_new=prior.new, runif1=runif1,
                    betas=betas, phi=phi,ntot=ntot,
                    ncomm=ncomm,nparam=nparam, nloc=nloc, nspp=nspp)
  
  list(betas=betas,accept=betas!=betas.orig)
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