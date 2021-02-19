#' @title FAB inference for generalized linear models
#' 
#' @description asymptotic FAB p-values and confidence intervals for parameters
#' in generalized linear regression models  
#' 
#' @param cformula formula for the control variables
#' @param FABvars  matrix of regressors for which to make FAB p-values and CIs
#' @param lformula formula for the linking model (just specify right-hand side)
#' @param alpha error rate for CIs (1-alpha CIs will be constructed) 
#' @param silent show progress (TRUE) or not (FALSE) 
#' @param ... additional arguments to be passed to \code{glm} 
#' 
#' @return an object of the class \code{glmFAB} which inherits from \code{glm}
#' 
#' @author Peter Hoff 
#'
#' @examples
#' 
#' # n observations, p FAB variables, q=2 control variables 
#' 
#' n<-100 ; p<-25 
#'
#' # X is design matrix for params of interest
#' # beta is vector of true parameter values 
#' # v a variable in the linking model - used to share info across betas
#' 
#' v<-rnorm(p) ; beta<-(2 - 2*v + rnorm(p))/3 ; X<-matrix(rnorm(n*p),n,p)/8
#' 
#' # control coefficients and variables  
#' alpha1<-.5 ; alpha2<- -.5
#' w1<-rnorm(n)/8
#' w2<-rnorm(n)/8
#' 
#' # simulate data 
#' lp<-1 + alpha1*w1 + alpha2*w2 + X%*%beta 
#' y<-rpois(n,exp(lp))
#' 
#' # fit model
#' fit<-glmFAB(y~w1+w2,X,~v,family=poisson)
#' 
#' fit$FABpv
#' fit$FABci 
#' summary(fit) # look at p-value column 
#'
#' @export 
glmFAB<-function(cformula, FABvars, lformula=NULL, alpha=.05,silent=FALSE,...){

  ## sampling model
  y<-model.frame(cformula)[[1]]
  W<-model.matrix(cformula) 

  ## the way lm is called below requires unique names
  wnames<-table(colnames(W)) 
  badnames<-names(wnames)[which(wnames!=1) ] 
  idx<-which( is.element( colnames(W), badnames) ) 
  colnames(W)[idx]<-paste0("v",1:length(idx))  

  X<-FABvars  
  p<-ncol(X) ; q<-ncol(W) 
  if(is.null(colnames(X))){ colnames(X)<-paste0("fv",1:p) }
 
  ## fit sampling model 
  if(!silent){ cat("\n") ; cat("Fitting sampling model: - ") }
  fit<-glm(y~.+0,as.data.frame(cbind(W,X)),...)  
  if(!silent){cat("\r","Fitting sampling model: # ","\n")}

  betaMLE<-fit$coef[q+1:p] 
  Sigma<-summary(fit)$cov.scaled[q+1:p,q+1:p]  
  s2<-diag(Sigma) 
  ZSTAT<-betaMLE/sqrt(s2) 
 
  ## fit linking model for each testvar 
  if(is.null(lformula)){ lformula<-formula( ~ -1 + rep(1,p) ) } 
  V<-model.matrix(lformula)

  RC<-T2<-NULL
  pbup<-min(50,p) 
  for(j in 1:p){ 

    ## progress
    pbprog<-floor(pbup*j/p)  
    pbar<-paste0(rep(c("#","-"),times=c(pbprog,pbup-pbprog)),collapse="")
    if(!silent){cat("\r","Fitting linking models:",pbar)}

    G<-MASS::Null(Sigma[,j]) 
    fitFH<-mmleFHP(t(G)%*%betaMLE, t(G)%*%V, t(G)%*%Sigma%*%G  )
    RC<-rbind(RC,fitFH$beta) ; T2<-c(T2,fitFH$t2) 
  }
  if(!silent){cat("\n") }

  Ebeta<-apply(RC*V,1,sum)
  Vbeta<-T2

  ## FAB p-values 
  b<-2*sqrt(s2)*Ebeta/Vbeta
  pFAB<-1 - abs( pnorm(ZSTAT+b) - pnorm(-ZSTAT) )

  ## FAB CI 
  CI<-NULL
  for(j in 1:p){
    CI<-rbind(CI,fabzCI(betaMLE[j],Ebeta[j],Vbeta[j],s2[j],alpha)) 
  } 

  ## put it together 
  rownames(CI)<-names(Ebeta)<-names(Vbeta)<-names(pFAB)
  fit$FABvars<-q+1:p
  fit$FABpv<-pFAB
  fit$FABci<-CI
  fit$Ebeta<-Ebeta
  fit$VBeta<-Vbeta 
  fit$SE<-sqrt(s2) 
  class(fit)<-c("glmFAB","glm") 
  fit
}
