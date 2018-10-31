plot(llk[9000:ngibbs],type='l')

#find number of groups detected
boxplot(theta)
ncomm=10
ind=c(1:3,9,10)
theta1=theta[,ind]

boxplot(theta1)
plot(NA,NA,ylim=c(0,1),xlim=c(0,nloc))
for (i in 1:ncomm) lines(theta1[,i],col=i)

par(mfrow=c(2,1),mar=rep(1,4))
plot(NA,NA,ylim=c(0,1),xlim=c(0,nloc))
for (i in 1:5) lines(theta1[,i],col=i)
plot(NA,NA,ylim=c(0,1),xlim=c(0,nloc))
for (i in 6:10) lines(theta1[,i],col=i)

#---------------------------------------------
#COMPARE LIKELIHOOD 
#it is complicated to compare log-prior because ncol(theta.true)=5 and not 10 and therefore sum(dbeta(vmat)) will be different

#true likelihood
setwd('U:\\presence absence model\\simulations') 
ncommun=5
nome=paste(c('fake data ','theta ','phi '),ncommun,'.csv',sep='')
theta.true=data.matrix(read.csv(nome[2],as.is=T))
phi.true=data.matrix(read.csv(nome[3],as.is=T))

prob=theta.true%*%phi.true
prob[prob>hi]=hi; prob[prob<lo]=lo
prob1=prob[loc.id,]
true.logl=sum(dat1*log(prob1)+(1-dat1)*log(1-prob1))
plot(llk,type='l',ylim=range(c(llk,true.logl)))
abline(h=true.logl,col='red')

