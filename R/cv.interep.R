#' This function does k-fold cross-validation for interep and returns the optimal value of lambda.
#' @importFrom stats gaussian
#' @importFrom MASS ginv
#' @param e matrix of environment factors.
#' @param z matrix of omics factors. In the case study, the omics measurements are lipidomics data.
#' @param y the longitudinal response.
#' @param response type of the longitudinal response, the default is continuous.
#' @param initiation the method for iniating the coefficient vector.  The default is lasso.
#' @param alpha.i the elastic-net mixing parameter.  The program adopts the elastic-net to choose initial values of the coefficient vector.
#' alpha.i is the elastic-net mixing parameter, with 0 \eqn{\le} alpha.i \eqn{\le} 1.  alpha.i=1 is the lasso penalty, and alpha.i=0 is the ridge penalty.
#' The default is 1.  If the user chooses a method other than elastic-net to initialize coefficients, alpha.i will be ignored.
#' @param lambda1 a user-supplied sequence of \eqn{\lambda_{1}} values, which serves as a tuning parameter for individual predictors.
#' @param lambda2 a user-supplied sequence of \eqn{\lambda_{2}} values, which serves as a tuning parameter for interactions.
#' @param nfolds the number of folds for cross-validation.
#' @param maxits the maximum number of iterations that is used in the estimation algorithm.  The default value is 30.
#' @param corre the working correlation structure that is used in the estimation algorithm. interep provides three choices for the
#' working correlation structure: "AR-1", "independece" and "exchangeable".
#' @details
#' When dealing with predictors with both main effects and interactions, this function returns two optimal tuning parameters,
#' \eqn{\lambda_{1}} and \eqn{\lambda_{2}}; when there are only main effects in the predictors, this function returns \eqn{\lambda_{1}},
#' which is the optimal tuning parameter for individual predictors containing main effects.
#' @return an object of class "cv.interep" is returned, which is a list with components:
#' \item{lam1}{the optimal \eqn{\lambda_{1}}.}
#' \item{lam2}{the optimal \eqn{\lambda_{2}}.}
#' @references
#' Zhou, F., Wang, W., Jiang, Y. and Wu, C. (2018+). Variable selection for interactions in longitudinal lipidomics studies.
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


cv.interep <- function(e, z, y, response="continuous", initiation=NULL, alpha.i=1,
                       lambda1, lambda2, nfolds, maxits=30, corre){
  n=dim(y)[1]
  k=dim(y)[2]
  q=dim(e)[2]
  pindex=c(1:(q+1))
  p1=dim(z)[2]
  x=cbind(e,z)
  for (i in 1:p1) {
    for (j in 1:q) {
      x=cbind(x,e[,j]*z[,i])
    }
  }

  #==========================================find initial values for beta using glmnet==========================#

  x1=cbind(data.frame(rep(1,n)),x)
  x1=data.matrix(x1)
  lasso.cv <- glmnet::cv.glmnet(x1,y[,3],alpha=alpha.i,nfolds=5)
  alpha <- lasso.cv$lambda.min/10  # lambda in the notes
  # alpha <- 0.5
  lasso.fit <- glmnet::glmnet(x1,y[,3],family="gaussian",alpha=alpha.i,nlambda=100)
  beta.new <- as.vector(stats::predict(lasso.fit, s=alpha, type="coefficients"))[-1]

  #==========================================reformat the data===============================#
  data=reformat(k,y,x)
  y=data$y
  x=data$x
  id=data$id

  p=dim(x)[2]

  l=length(lambda1)


  ll=length(lambda2)


  k=rep(k,n)
  aindex=cumsum(k)
  index=c(0,aindex[-length(aindex)])
  n.train=(nfolds-1)/nfolds*n
  eps=0.001

  if(response=="continuous"){
    family = gaussian(link = "identity")
  }
  else if (response==NULL){
    family = gaussian(link = "identity")
  }
  #pay attention to this part
  lam1.min <- -1
  lam2.min <- -1
  cv.min <- Inf
  cv.vect <- NULL
  #=========================================GEE MCP group MCP======================================#
  if(q>1){
    for (l1 in 1:l) {

      lam1=lambda1[l1]
      for (l2 in 1:ll) {

        #lam2=lam1*t12
        lam2=lambda2[l2]

        sse=0
        #Perform nfolds cross validation
        for(cv in 1:nfolds){
          #Segement your data by fold using the which() function
          testIndexes <- ((cv-1)*k[1]*(n/nfolds)+1):(cv*k[1]*(n/nfolds))
          x.test <- x[testIndexes, ]
          y.test <- y[testIndexes]
          x.train <- x[-testIndexes, ]
          y.train <- y[-testIndexes]
          converge=F
          iter=0

          while ((!converge) & (iter < maxits)) {
            beta = beta.new
            Rhat=CorrR(y.train,x.train,beta,n,k,corre)
            Score=ScoreU(n,k,y.train,x.train,p,response,beta,Rhat)
            U=Score$U
            qU=Score$qU

            beta.mcp=beta[1:(p1+q+1)]
            x.mcp=x.train[,1:(p1+q+1)]

            E.mcp=rep(0,(p1+q+1))

            for (j in 1:(p1+q+1)) {
              sub=j
              x.sub=x.mcp[,sub]
              beta0=beta.mcp[sub]
              kj=t(x.sub)%*%x.sub/n
              norm = sqrt(mean((x.sub*beta0)^2))
              E.mcp[j]=dmcp(abs(as.vector(beta0)),lam1)/(abs(as.vector(norm))+eps)
            }

            x.gmcp=x.train[,(p1+q+2):p]
            beta.gmcp=beta[c((p1+q+2):p)]
            for (j in 1:p1) {
              sub=((j-1)*q+1):(j*q)
              x.sub=x.gmcp[,sub]
              beta0=beta.gmcp[sub]
              norm = sqrt(mean((x.sub%*%beta0)^2))
              E.mcp=c(E.mcp,dmcp(abs(as.vector(beta0)),lam2)/(abs(as.vector(norm))+eps))
            }

            E1.mcp=diag(E.mcp)
            E1.mcp[,pindex]<-0
            E.mcp<-E1.mcp

            E<-E1.mcp
            mat=qU + n*E
            mat[abs(mat)<0.000001]=0
            beta.new = beta + ginv(mat)%*%(U - n*E%*%beta)

            diff=mean(abs(beta.new-beta))
            converge = (diff < 1e-3)
            iter = iter+1
            #cat("iter",iter,"diff",diff,"\n")
          }
          eta=x.test%*%beta.new
          mu=family$linkinv(eta)
          ##family$dev.resids gives the square of the residuals
          sse <- sse+sum((family$dev.resids(y.test,mu,wt=1)))
        }
        cv.vect<-c(cv.vect, sse)

        if(sse<cv.min) {
          lam1.min<-lam1
          lam2.min <- lam2
          cv.min<-sse
        }
      }
    }
    return(list("lam1"=lam1,"lam2"=lam2))
  }
  #=========================================GEE MCP ======================================#
  else if(q==1){
    for (l1 in 1:l) {
      #lam1=0.85
      lam1=lambda1[l1]

      sse=0
      #Perform nfolds cross validation
      for(cv in 1:nfolds){
        #Segement your data by fold using the which() function
        testIndexes <- ((cv-1)*k[1]*(n/nfolds)+1):(cv*k[1]*(n/nfolds))
        x.test <- x[testIndexes, ]
        y.test <- y[testIndexes]
        x.train <- x[-testIndexes, ]
        y.train <- y[-testIndexes]
        converge=F
        iter=0

        while ((!converge) & (iter < maxits)) {
          beta = beta.new
          Rhat=CorrR(y.train,x.train,beta,n,k,corre)
          Score=ScoreU(n,k,y.train,x.train,p,response,beta,Rhat)
          U=Score$U
          qU=Score$qU

          beta.mcp=beta
          x.mcp=x.train

          E.mcp=rep(0,p)

          for (j in 1:p) {
            sub=j
            x.sub=x.mcp[,sub]
            beta0=beta.mcp[sub]
            kj=t(x.sub)%*%x.sub/n
            norm = sqrt(mean((x.sub*beta0)^2))
            E.mcp[j]=dmcp(abs(as.vector(beta0)),lam1)/(abs(as.vector(norm))+eps)
          }
          E1.mcp=diag(E.mcp)
          E1.mcp[,pindex]<-0
          E.mcp<-E1.mcp

          E<-E1.mcp
          mat=qU + n*E
          mat[abs(mat)<0.000001]=0
          beta.new = beta + ginv(mat)%*%(U - n*E%*%beta)

          diff=mean(abs(beta.new-beta))
          converge = (diff < 1e-3)
          iter = iter+1
          #cat("iter",iter,"diff",diff,"\n")
        }
        eta=x.test%*%beta.new
        mu=family$linkinv(eta)
        ##family$dev.resids gives the square of the residuals
        sse <- sse+sum((family$dev.resids(y.test,mu,wt=1)))
      }
      cv.vect<-c(cv.vect, sse)

      if(sse<cv.min) {
        lam1.min<-lam1
        lam2.min <- lam2
        cv.min<-sse
      }

    }
  }
  return(list("lam1"=lam1))
}
