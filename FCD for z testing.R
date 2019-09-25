loc=10
species=5

#assume that this are the results after having removed the i-th individual
nks=t(array.lsk[loc,,])
nlk=apply(array.lsk,c(1,3),sum)
media=1:ncomm

#calculate using full equations
prob=rep(0,ncomm)
for (ksel in 1:ncomm){
  nks.new=nks
  nlk.new=nlk
  
  nks.new[ksel,species]=nks.new[ksel,species]+1
  nlk.new[loc,ksel]=nlk.new[loc,ksel]+1
  
  tmp=sum(dpois(nlk.new[loc,],media,log=T))
  tmp1=0
  for (j in 1:ncomm){
    tmp1=tmp1+dmultinom(nks.new[j,],size=nlk.new[loc,j],prob=phi[j,],log=T)  
  }
  prob[ksel]=tmp1+tmp
}
prob=prob-max(prob)
prob=exp(prob)
prob=prob/sum(prob)
prob.full=prob

#calculate using subset of expressions
prob=rep(0,ncomm)
for (ksel in 1:ncomm){
  prob[ksel]=log(phi[ksel,species])+log(media[ksel])-log(nks[ksel,species]+1)
}
prob=prob-max(prob)
prob=exp(prob)
prob=prob/sum(prob)
prob.subset=prob
