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

// This function makes multiple multinomial draws
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

// This function makes a single multinomial
// [[Rcpp::export]]
int RmultinomSingle(double runif1, NumericVector prob, int ncommun) {
  double probcum = prob[0];
  int oo=0;
  
  while(runif1>probcum){
    oo=oo+1;
    probcum = probcum + prob[oo];
  }
  return oo;
}

// This function calculates ArrayP1a
// [[Rcpp::export]]

NumericVector GetArrayP1a(NumericMatrix lphi, NumericMatrix lmedia, int nloc, int nspp, int ncomm){
  arma::cube AY = arma::zeros<arma::cube>(nloc, nspp, ncomm);

  for (int l = 0; l < nloc; l++) {
    for (int s = 0; s < nspp; s++){
      for (int k = 0; k < ncomm; k++){
        AY(l,s,k) = lphi(k,s)+lmedia(l,k);
      }
    }
  }

  return(wrap(AY));
}

// This function samples z
// [[Rcpp::export]]

List SampleZ(IntegerMatrix nlk, IntegerMatrix yls, NumericVector ArrayP1,
             int nloc, int nspp, int ncomm, NumericVector ArrayLSK, int ntot){
  
  //initialize objects
  arma::cube ArrayLSKnew = arma::zeros<arma::cube>(nloc, nspp, ncomm);
  NumericVector lprob(ncomm);
  NumericVector prob(ncomm);
  NumericVector runif1=runif(ntot); //generate uniform random variables (I have to fix this *2)
  int ind;
  int oo=0;

  //convert array into arma::cube
  NumericVector vecArray(ArrayP1);
  arma::cube ArrayP1a(vecArray.begin(), nloc, nspp, ncomm, false);
  NumericVector vecArray1(ArrayLSK);
  arma::cube ArrayLSKa(vecArray1.begin(), nloc, nspp, ncomm, false);

  for (int l = 0; l < nloc; l++) {
    for (int s = 0; s < nspp; s++){
      
      //just proceed this if individuals for l and s exist
      if (yls(l,s)!=0){ 
        for (int k = 0; k < ncomm; k++){
          
          //just proceed if individuals for l, s, and k exist (logic problem right here)
          if (ArrayLSKa(l,s,k)!=0){ 
              for (int indiv = 0; indiv < ArrayLSKa(l,s,k); indiv++){ //loop for each individual
                //adjust number of individuals (i.e., calculate NlkStar)
                nlk(l,k)=nlk(l,k)-1;

                //calculate probabilities
                for (int k1 = 0; k1 < ncomm; k1++){
                  lprob[k1]=ArrayP1a(l,s,k1)-log(nlk(l,k1)+1);
                }
                lprob=lprob-max(lprob);
                lprob=exp(lprob);
                prob=lprob/sum(lprob);
                
                //sample z
                ind=RmultinomSingle(runif1[oo], prob, ncomm); 
                oo=oo+1;
                  
                //update matrices and arrays
                nlk(l,ind)=nlk(l,ind)+1;
                ArrayLSKnew(l,s,ind)=ArrayLSKnew(l,s,ind)+1;
              }
          }
        }
      }
    }
  }

  return Rcpp::List::create(Rcpp::Named("ArrayLSK") = ArrayLSKnew,
                            Rcpp::Named("nlk") = nlk);
}