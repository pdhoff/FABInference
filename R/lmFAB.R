#' @title FAB inference for linear models
#' 
#' @description FAB p-values and confidence intervals for parameters
#' in linear regression models  
#' 
#' @param cformula formua for the control variables
#' @param FABvars  matrix of regressors for which to make FAB p-values and CIs
#' @param lformula formula for the lining model (just specify right-hand side)
#' @param alpha error rate for CIs (1-alpha CIs will be constructed) 
#' @param rssSplit use some residual degrees of freedom to help fit linking model (TRUE/FALSE) 
#' @param silent show progress (TRUE) or not (FALSE) 
#' 
#' @examples
#' 
#' # n obervations, p FAB variables, q=2 control variables 
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
#' y<-rnorm(n,lp) 
#' 
#' # fit model
#' fit<-lmFAB(y~w1+w2,X,~v)
#' 
#' fit$FABpv
#' fit$FABci 
#' summary(fit) # look at p-value column 
#'
#' @export 
lmFAB<-function(cformula, FABvars, lformula=NULL, alpha=.05, rssSplit=TRUE, silent=FALSE ){

  ## sampling model
  mframe<-model.frame(cformula)
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
  fit<-lm(y~.+0,as.data.frame(cbind(W,X))) 
  betaOLS<-fit$coef[q+1:p]

  if(rssSplit==FALSE){ 
    df<-length(y) - (p+q) ; s2<-sum(fit$res^2)/df 
    df0<-0 ; ss0<-0
  }

  if(rssSplit!=FALSE){  
    if(0<rssSplit & rssSplit<1){ 
      df0<-floor(rssSplit*fit$df) 
    } else { df0=max(1,floor(fit$df/10))  }  
    rss<-rssSplit(fit,df0) 
    df<-fit$df-df0 
    s2<-rss[2]/df 
    ss0<-rss[1] 
  } 
  
  S<-solve(crossprod(cbind(W,X)))[q+1:p,q+1:p]
  SED<-sqrt(s2*diag(S))
  TSTAT<-betaOLS/SED
  if(!silent){cat("\r","Fitting sampling model: # ","\n")}

  ## fit linking model for each testvar 
  if(is.null(lformula)){ lformula<-formula( ~ -1 + rep(1,p) ) } 
  V<-model.matrix(lformula)

  RC<-T2<-S2<-NULL

  pbup<-min(50,p) 
  for(j in 1:p){ 

    ## progress
    pbprog<-floor(pbup*j/p)  
    pbar<-paste0(rep(c("#","-"),times=c(pbprog,pbup-pbprog)),collapse="")
    if(!silent){cat("\r","Fitting linking models:",pbar)} 

    G<-MASS::Null(S[,j])
    fitFH<-mmleFH(t(G)%*%betaOLS, t(G)%*%V, t(G)%*%S%*%G, ss0, df0)
    RC<-rbind(RC,fitFH$beta) ; T2<-c(T2,fitFH$t2) ; S2<-c(S2,fitFH$s2)  
  }
  if(!silent){cat("\n") }

  Ebeta<-apply(RC*V,1,sum)
  Vbeta<-T2
  SEI<-sqrt(S2*diag(S))

  ## FAB p-values 
  b<-2*SEI*Ebeta/Vbeta
  pFAB<-1 - abs( pt(TSTAT+b,df) - pt(-TSTAT,df) )

  ## FAB CI 
  CI<-NULL
  for(j in 1:p){
    CI<-rbind(CI,fabtzCI(betaOLS[j],SED[j],df,alpha,
                         list(mu=Ebeta[j],tau2=Vbeta[j],sigma2=SEI[j]^2)) )
  } 

  ## put it together 
  rownames(CI)<-names(Ebeta)<-names(Vbeta)<-names(pFAB)
  fit$FABvars<-q+1:p
  fit$FABpv<-pFAB
  fit$FABci<-CI
  fit$Ebeta<-Ebeta
  fit$VBeta<-Vbeta 
  fit$SEI<-SEI
  class(fit)<-c("lmFAB","lm") 
  fit
} 
 
 
 

