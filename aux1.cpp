#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <fstream>
using namespace Rcpp;

/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/

// This function makes multinomial draws
// [[Rcpp::export]]
IntegerVector rmultinom1(NumericVector runif1, NumericVector prob, int ncommun) {
  IntegerVector res(ncommun);
  double probcum = prob(0);
  
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

//' This function samples z's
// [[Rcpp::export]]
List samplez(NumericMatrix theta, NumericMatrix phi, IntegerMatrix y, int ncommun, int nloc, int nspp) {
  
  IntegerMatrix nlk(nloc,ncommun);
  IntegerMatrix nks(ncommun,nspp);
  
  NumericVector prob(ncommun);
  IntegerVector znew(ncommun);
  
  for(int i=0; i<nloc; i++){
    for (int j=0; j<nspp; j++){
      for (int k=0; k<ncommun; k++){
        prob(k)=theta(i,k)*phi(k,j);
      }
      prob=prob/sum(prob);

      //multinomial draw
      znew=rmultinom1(runif(y(i,j)),prob,ncommun);
      
      //add to results to export
      nlk(i,_)=nlk(i,_)+znew;
      nks(_,j)=nks(_,j)+znew;
    }
  }
  return List::create(Named("nlk") = nlk,
                      Named("nks") = nks);
}

//' This function converts vmat into theta
// [[Rcpp::export]]
NumericMatrix convertVtoTheta(NumericMatrix vmat,
                                NumericVector prod) {
  NumericMatrix res(vmat.nrow(),vmat.ncol());
  
  for(int j=0; j<vmat.ncol();j++){
    res(_,j)=vmat(_,j)*prod;    
    prod=prod*(1-vmat(_,j));
  }
  
  return (res);
}

//' This function calculates ngreater
// [[Rcpp::export]]
IntegerMatrix ngreater(IntegerMatrix nlk,int nloc, int ncommun){
  IntegerMatrix ngreater(nloc,ncommun);
  int oo=ncommun-1;
  IntegerVector tmp(nloc);

  while (oo>=0){
    tmp=tmp+nlk(_,oo);
    ngreater(_,oo)=tmp;
    oo=oo-1;
  }
  return ngreater;
}

//' This function samples from truncated normal distrib and stores the sufficient statistics
// [[Rcpp::export]]
NumericMatrix getw(IntegerMatrix ge, NumericMatrix psi, NumericMatrix pnorm1,IntegerMatrix nlk,
          int nloc, int ncomm){
  
  NumericMatrix soma(nloc,ncomm-1);

  for(int i=0; i<nloc;i++){
    for (int k=0; k<(ncomm-1); k++){
      NumericVector zneg=runif(ge(i,k+1),0.0,pnorm1(i,k));
      zneg=qnorm(zneg,psi(i,k),1.0);
      zneg=ifelse(zneg < -30,-30,zneg); //to avoid zneg being -Inf
      NumericVector zpos=runif(nlk(i,k),pnorm1(i,k),1.0);
      zpos=qnorm(zpos,psi(i,k),1.0);
      zpos=ifelse(zpos > 30,30,zpos); //to avoid zpos being Inf
      soma(i,k)=sum(zneg)+sum(zpos);
    }
  }
  
  return soma;
}

// // [[Rcpp::export]]
// NumericVector teste(NumericVector x){
//   
//   NumericVector zneg=ifelse(is_infinite(x)==1,-30.0,x);
// 
//   return zneg;
// }
