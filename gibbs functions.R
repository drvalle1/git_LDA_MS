#this function generates dirichlet random variables (1 one for each row of alpha)
rdirichlet1=function(alpha,ncomm,nspp){
  tmp=matrix(rgamma(n=ncomm*nspp,alpha,1),ncomm,nspp)
  soma=matrix(rowSums(tmp),ncomm,nspp)
  tmp/soma
}
#---------------------------------------------
sample.array=function(nloc,nspp,jump1,array.lsk,y,ncomm,media){
  accept.ls=matrix(0,nloc,nspp)
  
  for (i in 1:nloc){
    table.old=array.lsk[i,,]
    for (j in 1:nspp){
      table.new=table.old
      
      #propose new values for species j
      tmp=table.old[j,]+jump1[i,j]
      pprop.old=tmp/sum(tmp)
      table.new[j,]=rmultinom(1,size=y[i,j],prob=pprop.old)
      
      #calculate proposal lprobabilities
      lprob.old.to.new=ldmultinom1(x=table.new[j,],size=y[i,j],prob=pprop.old)
      tmp=table.new[j,]+jump1[i,j]
      pprop.new=tmp/sum(tmp)
      lprob.new.to.old=ldmultinom1(x=table.old[j,],size=y[i,j],prob=pprop.new)

      #calculate ltarget distribution
      ltarget.old=ltarget.new=0
      nk.old=colSums(table.old)
      nk.new=colSums(table.new)
      for (k in 1:ncomm){
        ltarget.old=ltarget.old+ldmultinom1(x=table.old[,k],size=nk.old[k],prob=phi[k,])
        ltarget.new=ltarget.new+ldmultinom1(x=table.new[,k],size=nk.new[k],prob=phi[k,])
      }

      #calculate lpriors
      lprior.old=sum(dpois(nk.old,media[i,],log=T))
      lprior.new=sum(dpois(nk.new,media[i,],log=T))
      
      #MH 
      p0=ltarget.old+lprior.old+lprob.old.to.new
      p1=ltarget.new+lprior.new+lprob.new.to.old
      k=acceptMH(p0=p0,p1=p1,x0=1,x1=2,T)
      if (k$x==2) {
        table.old=table.new
        accept.ls[i,j]=1
      }
    }
    array.lsk[i,,]=table.old
  }
  list(array.lsk=array.lsk,accept.ls=accept.ls)
}
#-------------------------------------------
ldmultinom1=function(size,x,prob){
  lgamma(size + 1) + sum(x * log(prob) - lgamma(x + 1))
}
sample.lambdas=function(lambda.a,lambda.b,nlk,ncomm,nloc,xmat,betas){
  nk=apply(nlk,2,sum)
  media=exp(xmat%*%betas)
  rgamma(ncomm,lambda.a+nk,lambda.b+apply(media,2,sum))
}
#----------------------------------------------
sample.betas=function(lambda.a,lambda.b,nlk,xmat,betas,ncomm,nparam,jump1){
  betas.orig=betas.old=betas.prop=betas
  betas.prop[]=rnorm(nparam*ncomm,mean=betas.old,sd=jump1)
  for (i in 1:nparam){
    betas.new=betas.old
    betas.new[i,]=betas.prop[i,]
    
    lmedia.old=xmat%*%betas.old
    lmedia.new=xmat%*%betas.new
    p1.old=apply(nlk*lmedia.old,2,sum)
    p1.new=apply(nlk*lmedia.new,2,sum)
    
    p2=lambda.a+apply(nlk,2,sum)
    p3.old=log(lambda.b+apply(exp(lmedia.old),2,sum))
    p3.new=log(lambda.b+apply(exp(lmedia.new),2,sum))
    
    p4.old=(1/2)*(betas.old[i,]^2)
    p4.new=(1/2)*(betas.new[i,]^2)
    
    pold=p1.old-p2*p3.old-p4.old
    pnew=p1.new-p2*p3.new-p4.new
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