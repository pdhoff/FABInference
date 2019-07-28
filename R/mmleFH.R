#' @title Marginal MLEs for the Fay-Herriot model
#' 
#' @description Marginal MLEs for the Fay-Herriot random effects model where 
#' the covariance matrix for the sampling model is known to scale. 
#' 
#' @param y direct data follwing normal model \eqn{y\sim N(\theta,V\sigma^2)} 
#' @param X linking model predictors \eqn{ \theta\sim N(X\beta,\tau^2 I)} 
#' @param V covariance matrix to scale
#' @param pSS prior sum of squares for estimate of \eqn{\sigma^2} 
#' @param pdf prior degrees of freedom for estimate of \eqn{\sigma^2} 
#' 
#' @export
mmleFH<-function(y,X,V,pSS=0,pdf=0){

  ## mml estimation under the model 
  ## $y \sim N(\theta,V \sigma^2)$ 
  ## $\theta \sim N( X\beta,\tau^2 I)$
  ## where $V$, $X$ are known. 
  ## Additional stability can be added by 
  ## including some information on \sigma^2 
  ## so that pSS/pdf \approx \sigma^2. 

  ## Also, could stabilize with pSS=pdf*sum(y^2)/sum(diag(V)), 
  ## which is centering around the null model of theta=0. 

  eV<-eigen(V)
  E<-eV$vec
  L<-eV$val

  obj<-function(lt2s2){
    t2<-exp(lt2s2[1])
    s2<-exp(lt2s2[2])
    G<-1/sqrt((t2+L*s2))

    yd<-G*crossprod(E,y)
    Xd<-G*crossprod(E,X) 

    RSS<-sum(lm(yd~ -1+Xd)$res^2)
    ldet<-2*sum(log(G))
    RSS - ldet  + pdf*log(s2) + pSS/s2
  }

  t2s2<-exp(optim(c(0,0),obj)$par) ; t2<-t2s2[1] ; s2<-t2s2[2]
  G<-1/sqrt((t2+L*s2))
  yd<-G*crossprod(E,y)
  Xd<-G*crossprod(E,X)
  fit<-lm(yd ~ -1+Xd)
  list(beta=fit$coef,t2=t2,s2=s2)
}


