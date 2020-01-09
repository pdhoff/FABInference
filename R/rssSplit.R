#' @title Residual sum of squares split
#' 
#' @description Split residual sum of squares from normal linear regression
#' 
#' @param fit lm object 
#' @param df0 degrees of freedom for the smaller of the two residual sums of squares 
#' @param seed random seed for constructing the basis vectors of the split
#' 
#' @return a two-dimensional vector of independent sums of squares
#' 
#' @author Peter Hoff 
#' 
#' @examples 
#' n<-30 ; p<-6 ; sigma2<-1.5 
#' X<-matrix(rnorm(n*p),n,p) 
#' y<-X%*%rnorm(6) + sqrt(sigma2)*rnorm(n) 
#' ss<-rssSplit(lm(y~ -1+X)) 
#' df<-as.numeric( substring(names(ss),first=3)) 
#' ss/df
#' 
#' @export
rssSplit<-function(fit,df0=max(1,floor(fit$df/10)),seed=-71407){

  e<-fit$res
  U<-qr.Q(fit$qr)
  n<-nrow(U)
  set.seed(seed)
  Z<-matrix(rnorm(n*df0),n,df0)
  V<-svd(Z-U%*%crossprod(U,Z))$u

  ss0<-sum( (t(V)%*%e)^2 )
  ss1<-sum(e^2) - ss0
  ss<-c(ss0,ss1)
  names(ss)<-paste0("df",c(df0,n-ncol(U)-df0) )
  ss
}




