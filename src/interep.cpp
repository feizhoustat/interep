#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include"reformat.h"
#include"CorrR.h"
#include"ScoreU.h"
#include"dmcp.h"
#include<Rmath.h>
#include<iostream>
#include<stdio.h>
#include <stdbool.h>
using namespace Rcpp;
using namespace arma;
using namespace R;

//' fit a generalized estimating equation with given lambda
//'
//' This function makes predictions for generalized estimating equation with a given value of lambda.
//' Typical usage is to have the cv.interep function compute the optimal lambda, then provide it to
//' the interep function.
//' @param e matrix of environment factors.
//' @param z matrix of omics factors. In the case study, the omics measurements are lipidomics data.
//' @param y0 the longitudinal response.
//' @param beta the intial value for the coefficient vector.
//' @param lam1 the tuning parameter lambda1 for individual predictors.
//' @param lam2 the tuning parameter lambda2 for interactions.
//' @param corre the working correlation structure that is used in estimation algorithm. interep provides three choices
//' for the working correlation structure: "a" as "AR-1", "i" as "independence" and "e" as "exchangeable".
//' @param maxits the maximum number of iterations that is used in the estimation algorithm.
//' @details
//' When dealing with predictors with both main effects and interactions, this function requires two optimal tuning parameters,
//' \eqn{\lambda_{1}} and \eqn{\lambda_{2}}; when there are only main effects in the predictors, this function only requires \eqn{\lambda_{1}},
//' @return
//' \item{coef}{the coefficient vector.}
//' @references
//' Zhou, F., Ren, J., Li X., Wang, W., Jiang, Y. and Wu, C. (2018+). Variable selection for interactions in longitudinal lipidomics studies.
//'
//' Wu, C., Zhong, P. & Cui, Y. (2018). Additive varying-coefficient model for nonlinear gene-environment interactions.
//' \href{https://doi.org/10.1515/sagmb-2017-0008}{\emph{Statistical Applications in Genetics and Molecular Biology}, 17(2)}
//'
//' Wu, C., Jiang, Y., Ren, J., Cui, Y. and Ma, S. (2018). Dissecting gene-environment interactions: a penalized robust approach accounting for hierarchical structures.
//' \href{https://doi.org/10.1002/sim.7518}{\emph{Statistics in Medicine}, 37:437–456}
//'
//' Wu, C., Shi, X., Cui, Y. and Ma, S. (2015) A penalized robust semiparametric approach for gene-environment interactions.
//' \href{https://doi.org/10.1002/sim.6609}{\emph{Statistics in Medicine}, 34 (30): 4016–4030}
//'
//' Wu, C., Cui, Y. and Ma, S. (2014) Integrative analysis of gene-environment interactions under a multi-response partially linear varying coefficient model.
//' \href{https://doi.org/10.1002/sim.6287}{\emph{Statistics in Medicine}, 33 (28): 4988–4498}
//'
//' Wu, C. and Cui Y. (2013) A novel method for identifying nonlinear gene-environment interactions in case-control association studies.
//' \href{https://doi.org/10.1007/s00439-013-1350-z}{\emph{Human Genetics}, 132 (12): 1413–1425}
//'
//' @examples
//' data("dat")
//' e=dat$e
//' z=dat$z
//' y=dat$y
//' beta=dat$beta
//' lam1=dat$lam1
//' lam2=dat$lam2
//' index=dat$index
//' b = interep(e, z, y, beta, lam1, lam2, corre='e',maxits=30)
//' b[abs(b)<0.05]=0
//' pos = which(b != 0)
//' tp = length(intersect(index, pos))
//' fp = length(pos) - tp
//' list(tp=tp, fp=fp)
//'
//' @export
// [[Rcpp::export()]]
arma::vec interep(arma::mat& e, arma::mat& z, arma::mat& y0, arma::vec& beta, double lam1, double lam2, char corre, int maxits){
  int n=y0.n_rows,
      k0=y0.n_cols,
      q=e.n_cols,
      p1=z.n_cols;
  arma::mat x0(n,(q+p1*(q+1))),
      x00(n,(q+p1*(q+1)));
  arma::vec betanew(beta.size()),
      a(n);

  for(int i=0; i<n; i++){
    for(int j=0; j<q; j++){
      x0(i,j)=e(i,j);
    }
  }

  for(int i=0; i<n; i++){
    for(int j=0; j<p1; j++){
      x0(i,(q+j))=z(i,j);
    }
  }

  for(int i=0; i<n; i++){
    for(int j=0; j<p1; j++){
      for(int jj=0; jj<q; jj++){
        x0(i,(q+p1+j*q+jj))=e(i,jj)*z(i,j);
      }
    }
  }

  //for(int i=0; i<(x0.n_cols); i++){
  for(int i=0; i<(q+p1*(q+1)); i++){
    for(int j=0; j<n; j++){
      a(j)=x0(j,i);
    }
    for(int j=0; j<n; j++){
      x00(j,i)=(x0(j,i)-mean(a))/stddev(a);
    }
  }
  //reformat the data
  Rcpp::List data=reformat(y0,x00);
  arma::vec y=data["y"];
  arma::mat x=data["x"];

  int p=x.n_cols;
  arma::vec k(n);
  k.fill(k0);

 //std::cout << k << std::endl;
  float diff,
        eps=0.001;
  bool converge=false;
  int iter=0;

  //interep
  arma::vec beta1(p1+q+1),
      beta2(3*p1),
      E1(p,arma::fill::zeros),
      E0(q),
      beta0(1),
      beta00(q);
  arma::mat x1(x.n_rows,(p1+q+1)),
      x2(x.n_rows,(3*p1)),
      xsub2(x.n_rows,q),
      E(p,p),
      matr(p,p);

  for(int i=0; i<(p1+q+1); i++){
    beta1(i)=beta(i);
  }

  for(int i=0; i<3*p1; i++){
    beta2(i)=beta(p1+q+1+i);
  }

  for(int i=0; i<n; i++){
    for(int j=0; j<(p1+q+1); j++){
      x1(i,j)=x(i,j);
    }
  }

  for(int i=0; i<n; i++){
    for(int j=0; j<3*p1; j++){
      x2(i,j)=x(i,(p1+q+1+j));
    }
  }


  while ((!converge) & (iter < maxits)) {
    arma::cube Rhat=CorrR(y,x,beta,n,k,corre);
    Rcpp::List Score=ScoreU(n,k,y,x,p,beta,Rhat);
    arma::vec U=Score["U"];
    arma::mat qU=Score["qU"];
    //std::cout << U << endl;
    for (int j=(q+1); j<(p1+q+1); j++) {
      int sub=j;
      arma::vec xsub1(x1.n_rows),
      kk(x1.n_rows);
      beta0(0)=beta1(sub);
      for(int i=0; i<(n*k0); i++){
      //for(int i=0; i<(x1.n_rows); i++){
        xsub1(i)=x1(i,sub);
        kk(i)=pow(xsub1(i)*beta0(0),2);
      }

      float norm = sqrt(mean(kk));
      E1(j)=dmcp(beta0,lam1,3)(0)/(std::abs(norm)+eps);
    }

    for (int j=0; j<p1; j++) {
      for(int jj=j*q; jj<(j+1)*q;jj++){
        for(int i=0; i<q; i++){
          beta00(i)=beta2(jj);
          for(int ii=0; ii<(n*k0); ii++){
          //for(int ii=0; ii<(x2.n_rows); ii++){
            xsub2(ii,i)=x2(ii,jj);
          }
        }
      }

      arma::vec eta=xsub2*beta00;
      double norm = sqrt(mean(eta%eta));
      E0=dmcp(beta00,lam2,3)/(std::abs(norm)+eps);
      for(int jj=j*q; jj<(j+1)*q;jj++){
        for(int i=0; i<q; i++){
          E1(1+q+p1+jj)=E0(i);
        }
      }
    }

    for(int i=0; i<p; i++){
      E(i,i)=n*E1(i);
    }

    for(int i=0; i<p; i++){
      for(int j=0; j<p; j++){
        matr(i,j)=qU(i,j)+E(i,j);
      }
    }

    betanew = beta + arma::inv(matr)*(U - E*beta);

    diff=mean(abs(betanew-beta));
    if(diff < 0.001){
      converge = true;
    }
    iter++;
    for(int i=0; i<p; i++){
      beta(i)=betanew(i);
    }
  }
  //std::cout <<iter << std::endl;
  //std::cout << diff << std::endl;

  return(betanew);
}






