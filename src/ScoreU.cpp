#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include"ScoreU.h"
#include<Rmath.h>
#include<iostream>
#include<stdio.h>

using namespace Rcpp;
using namespace R;

// [[Rcpp::export()]]
Rcpp::List ScoreU(int n, arma::vec& k, arma::vec& y, arma::mat& x, int p, arma::vec& beta, arma::cube& Rhat){
  arma::vec aindex=cumsum(k),
      index(aindex.size(),arma::fill::zeros),
      U(p,arma::fill::zeros);
  int len=aindex.size()-1;

  for(int i=0; i<len; i++){
    index(i+1)=aindex(i);
  }

  arma::mat qU=arma::mat(p,p,arma::fill::zeros);
  for(int i=0; i<n; i++){
    arma::mat A=arma::mat(k(i),k(i),arma::fill::eye);
    arma::vec res(k(i),arma::fill::zeros);
    arma::mat D=arma::mat(k(i),x.n_cols);
    for(int j=0; j<k(i); j++){
      for(int jj=0; jj<p; jj++){
        D(j,jj)=x((j+index(i)),jj);
      }
    }
    arma::vec mu=D*beta;
    for(int t=0; t<k(i); t++){
      res(t)=y(t+index(i))-mu(t);
    }
    arma::mat V=Rhat.slice(i);
    U+=D.t()*inv(V)*res;
    qU+=D.t()*inv(V)*D;
  }

  return Rcpp::List::create(Rcpp::Named("U") = U,
                            Rcpp::Named("qU") = qU);

}


