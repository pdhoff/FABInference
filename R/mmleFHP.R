#' @title Marginal MLEs for the Fay-Herriot model with known covariance
#' 
#' @description Marginal MLEs for the Fay-Herriot random effects model where 
#' the covariance matrix for the sampling model is known
#' 
#' @param y direct data follwing normal model \eqn{y\sim N(\theta,\Sigma)} 
#' @param X linking model predictors \eqn{ \theta\sim N(X\beta,\tau^2 I)} 
#' @param Sigma covariance matrix in sampling model
#' 
#' @export
mmleFHP<-function(y,X,Sigma){

  ## mml estimation under the model 
  ## $y \sim N(\theta,\Sigma)$ 
  ## $\theta \sim N( X\beta,\tau^2 I)$
  ## where $\Sigma$, $X$ are known. 

  eS<-eigen(Sigma)
  E<-eS$vec
  L<-eS$val

  obj<-function(t2){
    G<-1/sqrt(t2+L)

    yd<-G*crossprod(E,y)
    Xd<-G*crossprod(E,X) 

    RSS<-sum(lm(yd~ -1+Xd)$res^2)
    ldet<-2*sum(log(G))
    RSS - ldet  
  }

  t2<-optimize(obj,c(0,mean(y^2)))$min
  G<-1/sqrt(t2+L)
  yd<-G*crossprod(E,y)
  Xd<-G*crossprod(E,X)
  fit<-lm(yd ~ -1+Xd)
  list(beta=fit$coef,t2=t2)
}



