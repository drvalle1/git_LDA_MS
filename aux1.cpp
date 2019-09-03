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

// This function makes a multinomial draw for size>1
// [[Rcpp::export]]
IntegerVector rmultinom1(NumericVector prob, int size) {
  IntegerVector res(prob.length());
  double probcum = prob(0);
  NumericVector runif1=runif(size);
  NumericVector runif2 = Rcpp::clone(runif1);
  runif2.sort(false);
  int oo=0;
  
  for (int i = 0; i < runif2.length(); i++) {
    if (runif2(i)< probcum){
      res(oo)=res(oo)+1;
    } else {
      while (runif2(i)>probcum){
        oo=oo+1;
        probcum = probcum + prob(oo);  
      }
      res(oo)=res(oo)+1;
    }
  }
  return res;
}

// This function calculates the multinomial distribution
// [[Rcpp::export]]
double ldmultinom(IntegerVector x, int size, NumericVector prob) {
  NumericVector x1(x.length());
  
  //convert from integer to numeric
  for (int i = 0; i < x.length(); i++){
    x1[i]=(double) x[i];
  } 
  
  NumericVector res1=x1*log(prob);
  double res=lgamma(size+1)+sum(res1 - lgamma(x+1));
  return(res);
}


// this function calculates the loglikelih based on multinomial
// [[Rcpp::export]]
double LogLikMultin(int nloc,int ncomm, int nspp, NumericMatrix phi, NumericVector Arraylsk){
  double p2=0;
  //convert array into arma::cube
  NumericVector vecArray(Arraylsk);
  arma::cube ArrayLSK1(vecArray.begin(), nloc, nspp, ncomm, false);
  IntegerVector tmp(nspp);
  
  for (int l = 0; l < nloc; l++) {
    for (int k = 0; k < ncomm; k++){
      for (int s = 0; s < nspp; s++){
        tmp[s]=ArrayLSK1(l,s,k);
      }
      p2=p2+ldmultinom(tmp,sum(tmp),phi(k,_));
    }
  }
  return(p2);
}

// This function calculates log of Poisson probability
// [[Rcpp::export]]
double ldpois1(int x,double lambda){
  double res;
  res=x*log(lambda)-lambda-lgamma(x+1);
  return(res);
}

// This function implements the MH step
// [[Rcpp::export]]
int RcppAcceptMH(double lpOld,double lpNew, double runi){
  double a=exp(lpNew-lpOld);
  if (runi < a) return(2);
  else return(1);
}

// This function samples Array.lsk
// [[Rcpp::export]]

List SampleArray(NumericVector Arraylsk,int nloc, int nspp, int ncomm,
                          NumericMatrix jump1,IntegerMatrix y, 
                          NumericMatrix phi, NumericMatrix media,
                          NumericMatrix runif1){
  //convert array into arma::cube
  NumericVector vecArray(Arraylsk);
  arma::cube ArrayLSK1(vecArray.begin(), nloc, nspp, ncomm, false);
  
  //initialize stuff
  IntegerMatrix AcceptLS(nloc,nspp);
  IntegerMatrix TableOld(nspp,ncomm);
  IntegerMatrix TableNew(nspp,ncomm);
  NumericVector tmp(ncomm);
  NumericVector PpropNew(ncomm);
  NumericVector PpropOld(ncomm);
  double lprobOldtoNew;
  double lprobNewtoOld;
  double ltargetOld;
  double ltargetNew;
  double lpriorOld;
  double lpriorNew;
  double lpOld;
  double lpNew;
  int k;
  IntegerVector nkOld(ncomm);
  IntegerVector nkNew(ncomm);
  
  for (int l = 0; l < nloc; l++) {
    //fill-in TableOld based on ArrayLSK1
    for (int s = 0; s < nspp; s++){
      for (int k = 0; k < ncomm; k++){
        TableOld(s,k)=ArrayLSK1(l,s,k);
      }
    }
    
    for (int s = 0; s < nspp; s++){
      if (y(l,s)>0){
        TableNew=clone(TableOld);
        
        //propose new values for species s
        for (int k = 0; k < ncomm; k++){
          tmp[k]=TableOld(s,k)+jump1(l,s);  
        }
        PpropOld=tmp/sum(tmp);
        TableNew(s,_)=rmultinom1(PpropOld, y(l,s));
        
        //calculate proposal lprobabilities
        lprobOldtoNew=ldmultinom(TableNew(s,_),y(l,s),PpropOld);
        for (int k = 0; k < ncomm; k++){
          tmp[k]=TableNew(s,k)+jump1(l,s);  
        }
        PpropNew=tmp/sum(tmp);
        lprobNewtoOld=ldmultinom(TableOld(s,_),y(l,s),PpropNew);
        
        //calculate target and prior lprobabilities
        ltargetOld=0;
        ltargetNew=0;
        lpriorOld=0;
        lpriorNew=0;
        for (int k = 0; k < ncomm; k++){
          nkOld[k]=sum(TableOld(_,k));
          nkNew[k]=sum(TableNew(_,k));
          ltargetOld=ltargetOld+ldmultinom(TableOld(_,k),nkOld[k],phi(k,_));
          ltargetNew=ltargetNew+ldmultinom(TableNew(_,k),nkNew[k],phi(k,_));
          lpriorOld=lpriorOld+ldpois1(nkOld[k],media(l,k));
          lpriorNew=lpriorNew+ldpois1(nkNew[k],media(l,k));
        }
        
        //MH: accept or reject proposal
        lpOld=ltargetOld+lpriorOld+lprobOldtoNew;
        lpNew=ltargetNew+lpriorNew+lprobNewtoOld;
        k=RcppAcceptMH(lpOld,lpNew,runif1(l,s));
        
        if (k==2){
          TableOld=clone(TableNew);
          AcceptLS(l,s)=1;
        }
      }
    }
    
    //fill-in ArrayLSK1 based on TableOld
    for (int s = 0; s < nspp; s++){
      for (int k = 0; k < ncomm; k++){
        ArrayLSK1(l,s,k)=TableOld(s,k);
      }
    }
  }
  List L = List::create(Named("ArrayLSK") =ArrayLSK1, Named("AcceptLS") = AcceptLS);
  
  return(L);
}
