#' @useDynLib interep, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' k-folds cross-validation for interep
#'
#' This function does k-fold cross-validation for interep and returns the optimal value of lambda.
#' @param e matrix of environment factors.
#' @param z matrix of omics factors. In the case study, the omics measurements are lipidomics data.
#' @param y the longitudinal response.
#' @param beta the intial value for the coefficient vector.
#' @param lambda1 a user-supplied sequence of \eqn{\lambda_{1}} values, which serves as a tuning parameter for individual predictors.
#' @param lambda2 a user-supplied sequence of \eqn{\lambda_{2}} values, which serves as a tuning parameter for interactions.
#' @param nfolds the number of folds for cross-validation.
#' @param corre the working correlation structure that is used in the estimation algorithm. interep provides three choices for the
#' working correlation structure: "a" as AR-1", "i" as "independence" and "e" as "exchangeable".
#' @param maxits the maximum number of iterations that is used in the estimation algorithm.
#' @details
#' When dealing with predictors with both main effects and interactions, this function returns two optimal tuning parameters,
#' \eqn{\lambda_{1}} and \eqn{\lambda_{2}}; when there are only main effects in the predictors, this function returns \eqn{\lambda_{1}},
#' which is the optimal tuning parameter for individual predictors containing main effects.
#' @return an object of class "cv.interep" is returned, which is a list with components:
#' \item{lam1}{the optimal \eqn{\lambda_{1}}.}
#' \item{lam2}{the optimal \eqn{\lambda_{2}}.}
#' @references
#' Zhou, F., Ren, J., Li X., Wang, W., Jiang, Y. and Wu, C. (2018+). Variable selection for interactions in longitudinal lipidomics studies.
#'
#' Wu, C., Zhong, P. & Cui, Y. (2018). Additive varying-coefficient model for nonlinear gene-environment interactions.
#' \href{https://doi.org/10.1515/sagmb-2017-0008}{\emph{Statistical Applications in Genetics and Molecular Biology}, 17(2)}
#'
#' Wu, C., Jiang, Y., Ren, J., Cui, Y. and Ma, S. (2018). Dissecting gene-environment interactions: a penalized robust approach accounting for hierarchical structures.
#' \href{https://doi.org/10.1002/sim.7518}{\emph{Statistics in Medicine}, 37:437–456}
#'
#' Wu, C., Shi, X., Cui, Y. and Ma, S. (2015) A penalized robust semiparametric approach for gene-environment interactions.
#' \href{https://doi.org/10.1002/sim.6609}{\emph{Statistics in Medicine}, 34 (30): 4016–4030}
#'
#' Wu, C., Cui, Y. and Ma, S. (2014) Integrative analysis of gene-environment interactions under a multi-response partially linear varying coefficient model.
#' \href{https://doi.org/10.1002/sim.6287}{\emph{Statistics in Medicine}, 33 (28): 4988–4498}
#'
#' Wu, C. and Cui Y. (2013) A novel method for identifying nonlinear gene-environment interactions in case-control association studies.
#' \href{https://doi.org/10.1007/s00439-013-1350-z}{\emph{Human Genetics}, 132 (12): 1413–1425}
#'
#' @export



cv.interep <- function(e, z, y, beta, lambda1, lambda2, nfolds, corre, maxits){
  n=dim(y)[1]
  q=dim(e)[2]
  p1=dim(z)[2]
  len1=length(lambda1)
  len2=length(lambda2)

  folds=rep(1,n/nfolds)
  for (i in 2:nfolds) {
    folds=c(folds,rep(i,n/nfolds))
  }
  folds=sample(folds)
  pred=matrix(0,len1,len2)
  for (i in 1:len1) {
    lam1=lambda1[i]
    for (j in 1:len2){
      lam2=lambda2[j]
      mse=0
      for (cv in 1:nfolds) {
        #Segement your data by fold using the which() function
        testIndexes <- which(folds==cv,arr.ind=TRUE)
        e.test=e[testIndexes,]
        z.test=z[testIndexes,]
        y.test=y[testIndexes,]
        e.train=e[-testIndexes,]
        z.train=z[-testIndexes,]
        y.train=y[-testIndexes,]
        e.train=as.matrix(e.train)
        z.train=as.matrix(z.train)
        y.train=as.matrix(y.train)
        beta=interep(e.train, z.train, y.train, beta,lam1, lam2, corre, maxits)
        x.test=cbind(e.test,z.test)
        for (i1 in 1:p1) {
          for (j1 in 1:q) {
            x.test=cbind(x.test,e.test[,j1]*z.test[,i1])
          }
        }
        x.test=scale(x.test)
        data.test=reformat(y.test, x.test)
        x.test=data.test$x
        y.test=data.test$y
        mu=x.test%*%beta
        mse=mse+mean((y.test-mu)^2)
      }
      pred[i,j]=mse/nfolds
    }
  }
  lamb1=lambda1[which(pred==min(pred),arr.ind = TRUE)[1]]
  lamb2=lambda2[which(pred==min(pred),arr.ind = TRUE)[2]]
  return(list("lam1"=lamb1,"lam2"=lamb2))
}
