#this function generates dirichlet random variables (1 one for each row of alpha)
rdirichlet1=function(alpha,ncomm,nspp){
  tmp=matrix(rgamma(n=ncomm*nspp,alpha,1),ncomm,nspp)
  soma=matrix(rowSums(tmp),ncomm,nspp)
  tmp/soma
}
#---------------------------------------------
samplez.R=function(lphi, lmedia, array.lsk,nlk, ncomm,nloc, nspp,y){
  #calculate part 1 of the log-probability of z
  
  # teste=array(NA,dim=c(nloc,nspp,ncomm))
  # for (l in 1:nloc){
  #   for (s in 1:nspp){
  #     for (k in 1:ncomm){
  #       teste[l,s,k]=lphi[k,s]+lmedia[l,k]
  #     }
  #   }
  # }
  tmp = GetArrayP1a(lphi=lphi,lmedia=lmedia,nloc=nloc,nspp=nspp,ncomm=ncomm)
  array.p1 = array(tmp,dim=c(nloc,nspp,ncomm))
  # sum(array.p1!=teste)
  # identical(array.p1,teste)

  #sample z
  set.seed(1)
  teste=SampleZ(nlk=nlk, yls=y, ArrayP1=array.p1,
                nloc=nloc, nspp=nspp, ncomm=ncomm, ArrayLSK=array.lsk,ntot=sum(nlk))
  # sum(array.lsk!=teste$ArrayLSK) #things have changed!
  # identical(array.lsk,teste$ArrayLSK)
  #sum(teste$ArrayLSK); sum(y)
  
  # set.seed(1)
  # for (l in 1:nloc){
  #   nlk.sel=nlk[l,]
  # 
  #   for (s in 1:nspp){
  # 
  #     if (y[l,s]!=0){ #there is no point in sampling z if 0 individuals
  #       array.p1a=array.p1[l,s,]
  #       qtt=array.lsk[l,s,]
  #       indic=which(qtt!=0) #just work with those individuals that truly exist
  #       
  #       for (k in indic){
  #         qtt1=qtt[k]
  #         runif1=runif(qtt1)
  #         for (indiv in 1:qtt1){
  #           #adjust number of individuals (i.e., calculate n.star's)
  #           nlk.star=nlk.sel
  #           nlk.star[k]=nlk.star[k]-1
  # 
  #           #sample z
  #           lprob=array.p1a-log(nlk.star+1)
  #           lprob1=lprob-max(lprob)
  #           tmp=exp(lprob1)
  #           prob=tmp/sum(tmp)
  #           ind=RmultinomSingle1(runif1[indiv], prob, ncomm)
  # 
  #           #update matrices
  #           nlk.star[ind]=nlk.star[ind]+1
  #           nlk.sel=nlk.star
  # 
  #           #update arrays
  #           array.lsk[l,s,k]=array.lsk[l,s,k]-1
  #           array.lsk[l,s,ind]=array.lsk[l,s,ind]+1
  #         }
  #       }
  #     }
  #   }
  #   nlk[l,]=nlk.sel
  # }
  # sum(array.lsk!=teste$ArrayLSK) (cpp yields same results as R)
  # identical(array.lsk,teste$ArrayLSK)
  
  #basic test
  # sum(apply(teste$ArrayLSK,1:2,sum)!=y)
  # sum(apply(teste$ArrayLSK,c(1,3),sum)!=teste$nlk)

  list(nlk=teste$nlk,array.lsk=teste$ArrayLSK)        
}
#-------------------------------------------
sample.lambdas=function(lambda.a,lambda.b,nlk,ncomm,nloc){
  nk=apply(nlk,2,sum)
  rgamma(ncomm,lambda.a+nk,lambda.b+nloc)
}
#----------------------------------------------
RmultinomSingle1=function(runif1, prob, ncomm){
  probcum = prob[1];
  oo=1;
  
  while(runif1>probcum){
    oo=oo+1;
    probcum = probcum + prob[oo];
  }
  oo;
}
