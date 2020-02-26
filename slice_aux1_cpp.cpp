// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <fstream>
using namespace Rcpp;

/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/

// This function calculates the logtarget distribution
// [[Rcpp::export]]
double LogTarget(NumericVector LogMediaMiss, NumericVector param, NumericVector y, NumericMatrix xmat, 
                 int target, NumericVector var1) {
  NumericVector LogMediaComp=LogMediaMiss+xmat(_,target)*param[target];
  NumericVector media=exp(LogMediaComp);
  NumericVector LogPrior=(-1/(2*var1))*(param*param);
  double res=sum(y*LogMediaComp-media)+sum(LogPrior);
  return res;
}

//this function doubles the interval until we are outside the slice
// [[Rcpp::export]]
NumericVector doubling(NumericVector LogMediaMiss, double yslice, double w, NumericVector param,
                       NumericVector y, NumericMatrix xmat,int target, NumericVector var1){
  NumericVector ParamLo=clone(param);
  NumericVector ParamHi=clone(param);
  ParamLo[target]=ParamLo[target]-w*runif(1)[0];
  ParamHi[target]=ParamLo[target]+w;
  double ylo=LogTarget(LogMediaMiss,ParamLo,y,xmat,target,var1);
  double yhi=LogTarget(LogMediaMiss,ParamHi,y,xmat,target,var1);

  while(ylo>yslice){
    ParamLo[target]=ParamLo[target]-w;
    ylo=LogTarget(LogMediaMiss,ParamLo,y,xmat,target,var1);
  }
  while(yhi>yslice){
    ParamHi[target]=ParamHi[target]+w;
    yhi=LogTarget(LogMediaMiss,ParamHi,y,xmat,target,var1);
  }
  NumericVector res(2);
  res[0]=ParamLo[target];
  res[1]=ParamHi[target];
  return res;
}

//this function shrinks the slice if samples are outside the slice. If sample is inside the slice, accept this sample
// [[Rcpp::export]]
double SampleEachParam(NumericVector LogMediaMiss,NumericVector rango1,double yslice,NumericVector param,
                       NumericVector y,NumericMatrix xmat,int target,NumericVector var1) {
  double yfim=R_NegInf;
  double x=0;
  double DistLo;
  double DistHi;
  double diff1=rango1[1]-rango1[0];
  while ((yfim<yslice) & (diff1 > 0.00001)){
    x=rango1[0]+diff1*runif(1)[0];
    DistLo=abs(rango1[0]-x);
    DistHi=abs(rango1[1]-x);
    param[target]=x;
    yfim=LogTarget(LogMediaMiss,param,y,xmat,target,var1);
    if (yfim<yslice){
      if (DistLo<DistHi) rango1[0]=x;
      if (DistLo>DistHi) rango1[1]=x;
      diff1=rango1[1]-rango1[0];
    }
  }
  return x;
}

//this function performs matrix multiplication after excluding i-th column in xmat and i-th element in param
// [[Rcpp::export]]
NumericVector MatrixMultipExclude(NumericMatrix xmat, NumericVector param, int target, int nparam) {
  NumericVector res(xmat.nrow());
  for(int i=0;i<nparam;i++){
    if (i!=target){
      res=res+xmat(_,i)*param[i];
    }
  }
  return res;
}

//this function loops over all parameters for a given community, using a slice sampler for each parameter
// [[Rcpp::export]]
NumericVector SampleBetasEachComm(NumericVector param, NumericVector y, NumericMatrix xmat, 
                                  double w,int nparam, NumericVector var1) {
  NumericVector LogMediaMiss(xmat.nrow());
  double upper1;
  double yslice;
  NumericVector rango1(2);
  
  for(int i=0;i<nparam;i++){
    //define upper bound
    LogMediaMiss=MatrixMultipExclude(xmat,param,i,nparam);
    upper1=LogTarget(LogMediaMiss,param,y,xmat,i,var1);
    yslice=upper1-rexp(1)[0]; //method suggest by Neal 2003
    
    //define slice  
    rango1=doubling(LogMediaMiss,yslice,w,param,y,xmat,i,var1); //find range by doubling window
    
    //sample this particular parameter
    param[i]=SampleEachParam(LogMediaMiss,rango1,yslice,param,y,xmat,i,var1); //sample within the define range (rango1)
  }
  return param;
}

//this function samples all parameters for all communities
// [[Rcpp::export]]
NumericMatrix SampleBetas(NumericMatrix param, NumericMatrix y, NumericMatrix xmat, 
                          double w,int nparam,int ncomm, NumericVector var1) {
  NumericMatrix param1(nparam,ncomm);
  NumericVector res(nparam);
  for(int i=0;i<ncomm;i++){
    res=SampleBetasEachComm(param(_,i), y(_,i), xmat, w,nparam,var1);
    param1(_,i)=res;
  }
  return param1;
}

