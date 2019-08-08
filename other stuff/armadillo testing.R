library(RcppArmadillo)
library(inline)

if( require( RcppArmadillo ) ){
 	fx <- cxxfunction( signature(x = "matrix", n = "integer") , '
		Rcpp::NumericVector XX(x);
		Rcpp::IntegerVector dim(n);
		arma::cube AY = arma::zeros<arma::cube>(dim[0], dim[1], dim[2]);
		AY(0,0,0)=1;
		return(wrap(AY));
	', plugin = "RcppArmadillo" )
}
 
set.seed(10)
x = matrix(rnorm(10*15),10,15)
y = fx(x=x,n=c(25,25,10))
