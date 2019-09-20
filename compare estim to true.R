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

ordem=1:8#c(4,6,3,5,2,9,1,10) #4,6

#look at lambda
compare1(estim=res$lambda[ngibbs,ordem],true=lambda.true)

#look at nlk
tmp=matrix(res$nlk[ngibbs,],nloc,ncomm)
compare1(estim=tmp[,ordem],true=nlk.true)

#look at betas
k=matrix(res$betas[ngibbs,],nparam,ncomm)
compare1(estim=k[,ordem],true=betas.true)

#look at phi
tmp=matrix(res$phi[ngibbs,],ncomm,nspp)
tmp1=tmp[ordem,]
compare1(estim=tmp1,true=phi.true)

#look at convergence
# par(mfrow=c(4,3),mar=rep(1,4))
# for (i in 1:12){
#   plot(betas.out[,i],type='l')
# }
# 
# par(mfrow=c(3,2),mar=rep(1,4))
# for (i in 1:ncomm){
#   plot(lambda.out[,i],type='l')
# }
# 
# #look at jump
# z=jump1$array.lsk
# hist(z[y!=0]); table(z[y!=0])
# hist(z[y==0]); table(z[y==0])
# plot(jitter(y[y!=0]),jitter(jump1$array.lsk[y!=0]),ylim=c(0,20))