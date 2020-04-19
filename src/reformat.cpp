#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include"reformat.h"
#include<Rmath.h>
#include"reformat.h"
#include<iostream>
#include<stdio.h>

using namespace Rcpp;
using namespace arma;
using namespace R;

//' This function changes the format of the longitudinal data from wide format to long format
//' @param y the longitudinal response.
//' @param x a matrix of predictors, consisting of omics and environment factors, as well as their interactions. In
//' the case study, the omics measurements are lipidomics data.
//'
//' @examples
//' data("dat")
//' y=dat$y
//' x=dat$x
//' reformat(y,x)
//' @export
// [[Rcpp::export()]]
Rcpp::List reformat(arma::mat& y, arma::mat& x){
  int n=y.n_rows,
      k=y.n_cols,
      p1=x.n_cols;
  arma::vec y1(n*k);
  arma::mat x1=ones<mat>(n*k,1+p1);


  for(int i=0; i<n; i++){
    for(int t=0; t<k; t++){
      y1(i*k+t)=y(i,t);
      for(int j=0; j<p1; j++){
        x1(i*k+t,j+1)=x(i,j);
      }
    }
  }

    return Rcpp::List::create(Rcpp::Named("y") = y1,
                              Rcpp::Named("x") = x1);

}


