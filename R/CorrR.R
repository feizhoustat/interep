#' This function obtains the correlation matrix of longitudinal data
#' @importFrom stats sd
#' @param y the longitudinal response.
#' @param x a matrix of predictors, consisting of omics and environment factors, as well as their interactions. In
#' the case study, the omics measurements are lipidomics data.
#' @param beta the coefficient vector.
#' @param n the sample size.
#' @param k the number of repeated measurement.
#' @param corre the working correlation structure that is used in estimation algorithm. interep provides three choices
#' for the working correlation structure: "AR-1", "independence" and "exchangeable".
#' @return the correlation matrix.
#' @examples
#' data("dat")
#' y=dat$y
#' n=dim(y)[1]
#' x=dat$x
#' k=dat$k
#' data=reformat(k,y,x)
#' y=data$y
#' x=data$x
#' beta=dat$beta
#' corre="AR-1"
#' k=rep(dat$k,n)
#' CorrR(y,x,beta,n,k,corre)
#' @export

CorrR <- function(y,x,beta,n,k,corre){
  aindex=cumsum(k)
  index=c(0,aindex[-length(aindex)])
  mu=x%*%beta
  res=(as.vector(y)-mu)/sd(y)

  if(corre == "independence"){
    alfa_hat=0
  }
  else if(corre == "exchangeable"){
    res1=res[-c(rep(1,n)+c(0:(n-1))*k)]
    res2=res[-c((1:n)*k)]
    total=sum(res1*res2)
    alfa_hat=total/sum(k-1)

  }

  else if (corre == "AR-1") {
    res1=res[-c(rep(1,n)+c(0:(n-1))*k)]
    res2=res[-c((1:n)*k)]
    total=sum(res1*res2)
    alfa_hat=total/sum(k-1)

  }


  maxclsz=max(k)

  Rhat<-array(0,c(maxclsz,maxclsz,n))

  for (i in 1:n) {
    if (corre == "independence") {
      cor1 = diag(k[i])
    }
    else if (corre == "exchangeable") {
      cor1=matrix(alfa_hat,k[i],k[i])
      diag(cor1)=1
    }
    else if (corre == "AR-1"){
      cor1=matrix(0,k[i],k[i])
      for (t1 in 1:k[i]) {
        for (t2 in 1:k[i]) {
          cor1[t1,t2]<-alfa_hat^abs(t1-t2)
        }
      }
    }


    Rhat[1:k[i],1:k[i],i]<-cor1
  }

  return(Rhat)
}
