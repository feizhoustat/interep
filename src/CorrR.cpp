#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include"CorrR.h"
#include<Rmath.h>
#include<iostream>
#include<stdio.h>

using namespace Rcpp;
using namespace arma;
using namespace R;

//' This function obtains the correlation matrix of longitudinal data
//' @param y the longitudinal response.
//' @param x a matrix of predictors, consisting of omics and environment factors, as well as their interactions. In
//' the case study, the omics measurements are lipidomics data.
//' @param beta the coefficient vector.
//' @param n the sample size.
//' @param k a vector that contains the number of repeated measurements for each subject.
//' @param corre the working correlation structure that is used in estimation algorithm. interep provides three choices
//' for the working correlation structure: "a" as "AR-1", "i" as "independence" and "e" as "exchangeable".
//' @return the correlation matrix.
//' @examples
//' data("dat")
//' y=dat$y
//' n=dim(y)[1]
//' x=dat$x
//' data=reformat(y,x)
//' y=data$y
//' x=data$x
//' beta=dat$beta
//' corre='e'
//' k=rep(dat$k,n)
//' CorrR(y,x,beta,n,k,corre)
//' @export
// [[Rcpp::export()]]
arma::cube CorrR(arma::vec& y, arma::mat& x, arma::vec& beta, int n, arma::vec& k, char corre){
  arma::vec aindex=cumsum(k),
      index(aindex.size(),arma::fill::zeros);
  double alfa_hat=0;
  float sum5=0,
    sum6=0,
    sum7=0,
    sum8=0;
  int len=aindex.size()-1;
  for(int i=0; i<len; i++){
    index(i+1)=aindex(i);
  }
  arma::vec mu=x*beta;
 // vec res=(y-mu)/stddev(y);
 int len1=y.size();
 arma::vec res(len1);

 for(int i=0; i<len1; i++){
   res(i)=(y(i)-mu(i))/stddev(y);
 }

  if(corre == 'i'){
    alfa_hat=0;
  }
  else if(corre == 'e' || corre == 'a'){
 //else if(corre == 'e'){

    for (int i=0; i<n; i++){
      for (int j=0; j<k(i); j++){
        for (int jj=0; jj<k(i); jj++){
          if((j-jj)==1){
              sum7=res(j+index(i))*res(jj+index(i));
              sum5+=sum7;
          }
        }
      }
      sum8=(k(i)-1);
        sum6+=sum8;
    }
      alfa_hat=sum5/sum6;
  }


  int maxclsz=max(k);

  arma::cube Rhat(maxclsz,maxclsz,n);

  for (int i=0; i<n; i++) {
    arma::mat cor1 = arma::mat(k(i),k(i));
    if (corre == 'i') {
      cor1 = arma::mat(k(i),k(i),arma::fill::eye);
    }
    else if (corre == 'e') {
      cor1 = arma::mat(k(i),k(i),arma::fill::ones);
      for(int j=0; j<k(i); j++){
        for(int jj=0; jj<k(i); jj++){
          if(jj != j) {
            cor1(j,jj)=alfa_hat;
          }
        }
      }
    }
    else if (corre == 'a'){
      cor1 = arma::mat(k(i),k(i),arma::fill::ones);
      for(int j=0; j<k(i); j++){
        for(int jj=0; jj<k(i); jj++){
            cor1(j,jj)=pow(alfa_hat,std::abs(j-jj));
        }
      }
    }


    Rhat.slice(i)=cor1;
  }

  return(Rhat);

}


