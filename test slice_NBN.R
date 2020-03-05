rm(list=ls(all=TRUE))
library('Rcpp')
set.seed(10)

setwd('U:\\GIT_models\\git_LDA_MS')
sourceCpp('slice_NBN.cpp')

#generate fake data
nloc=100
nspp=50
media=2
Media=matrix(media,nloc,nspp)
n=10
tmp=rnbinom(nloc*nspp,mu=Media,size=n)
y=matrix(tmp,nloc,nspp)

#estimate NBN
ngibbs=1000
res=rep(NA,ngibbs)
NBN=1
w=10
#problem when i==9. Problem is with ShrinkAndSample
for (i in 1:ngibbs){
  print(i)
  NBN=SampleNBN(Media=Media,y=y,NBN=NBN,w=w)
  res[i]=NBN
}
plot(res,type='l')

#------------------------------
upper1=LogTargetNBN(Media=Media,y=y,NBN=NBN)
yslice=upper1-rexp(1); 
rango1=DoublingNBN(yslice=yslice, w=w,y=y,NBN=NBN,Media=Media)
rango2=ShrinkAndSample(Media=Media,rango1=rango1,yslice=yslice,y=y,NBN=NBN)
#------------------------------
upper1=LogTargetNBN(Media=Media,y=y,NBN=100)
upper2=LogTargetNBN(Media=Media,y=y,NBN=1)
upper1-upper2

upper1a=sum(dnbinom(y,mu=Media,size=100,log=T))
upper2a=sum(dnbinom(y,mu=Media,size=1,log=T))
upper1a-upper2a


sum(dnbinom(y,mu=Media,size=1000,log=T))
sum(dnbinom(y,mu=Media,size=10000,log=T))
LogTargetNBN(Media=Media,y=y,NBN=10000)