plot(res$llk[1:ngibbs],type='l')

compare1=function(estim,true){
  rango=range(c(true,estim))
  plot(true,estim,ylim=rango,xlim=rango)
  lines(rango,rango)  
}

k=res$lambda[ngibbs,]
plot(k,type='h')

k=res$betas[ngibbs,]
nparam=ncol(xmat)
k1=matrix(k,nparam,ncomm); round(k1,2)

ordem=c(5,3,4,2,1)
round(k1[,ordem],2)

#look at nlk
tmp=matrix(res$nlk[ngibbs,],nloc,ncomm); 
boxplot(tmp)
compare1(estim=jitter(tmp[,ordem]),true=jitter(nlk.true))

#look at betas
k=matrix(res$betas[ngibbs,],nparam,ncomm)
compare1(estim=k[,ordem],true=betas.true)

#look at phi
tmp=matrix(res$phi[ngibbs,],ncomm,nspp)
tmp1=tmp[ordem,]
compare1(estim=tmp1,true=phi.true)