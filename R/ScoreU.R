#' @importFrom stats gaussian
#' @importFrom MASS ginv
ScoreU <- function(n,k,y,x,p,response="continuous",beta,Rhat){
  #n=as.numeric(length(y)/k)
  aindex=cumsum(k)
  index=c(0,aindex[-length(aindex)])
  U=rep(0,p)
  qU=matrix(0,p,p)
  if(response=="continuous"){
    family = gaussian(link = "identity")
  }
  else if (response==NULL){
    family = gaussian(link = "identity")
  }
  for (i in 1:n){
    # D=rep(1,k) %*% t.default(x[i,])
    A=matrix(0,k[i],k[i])
    res=rep(0,k[i])
    D=x[c((1+index[i]):(k[i]+index[i])),]
    eta=D%*%beta
    mu=family$linkinv(eta)
    for (j in 1:k[i]) {
      # D=x[j+index[i],]
      # eta=D%*%beta
      # mu=family$linkinv(eta)
      A[j,j]=family$variance(y)[j+index[i]]
      #res[j]=y[i,j]-mu[j]
      res[j]=y[j+index[i]]-mu[j]
    }

    V=sqrt(A)%*%Rhat[1:k[i],1:k[i],i]%*%sqrt(A)
    U = U+t(D)%*%A%*%ginv(V)%*%res
    qU = qU + t(D)%*%A%*%ginv(V)%*%A%*%D
  }
  return(list("U"=U,"qU"=qU))
}
