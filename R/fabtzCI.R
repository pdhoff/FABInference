#' @title z-optimal FAB t-interval 
#' 
#' @description Computation of a 1-alpha FAB t-interval using 
#' z-optimal spending function
#' 
#' @param y a numeric scalar, a normally distributed statistic
#' @param s a numeric scalar, the standard error of y
#' @param dof positive integer, degrees of freedom for s
#' @param alpha the type I error rate, so 1-alpha is the coverage rate 
#' @param psi a list of parameters for the spending function, including 
#' \enumerate{
#' \item mu, the prior expectation of E[y]
#' \item tau2, the prior variance of E[y]
#' \item sigma2 the variance of y
#' }
#'
#' @examples
#' n<-10 
#' y<-rnorm(n) 
#' fabtzCI(mean(y),sqrt(var(y)/n),n-1)  
#' t.test(y)$conf.int 
#' @export
fabtzCI<-function(y,s,dof,alpha=.05,psi=list(mu=0,tau2=1e5,sigma2=1))
{
  mu<-psi$mu ; tau2<-psi$tau2 ; sigma2<-psi$sigma2
  if(tau2<=0)
  {
    thetaL<-min(mu,y+s*qt(alpha,dof))
    thetaU<-max(mu,y+s*qt(1-alpha,dof))
  }

  if(tau2>0)
  {

    root<-function(theta)
    {
      sfabz(theta,alpha=alpha,psi=psi) - pt( (y-theta)/s,dof )/alpha
    }
    a<-b<-y+s*qt(1-alpha,dof)
    #while(root(a)>0){ a<- a + s*qnorm(alpha)*.25 }
    #while(root(b)<0){ b<- b + s*qnorm(1-alpha)*.25 }  
    while(root(a)>0){ a<- a - s }
    while(root(b)<0){ b<- b + s }
    thetaU<-uniroot(root,c(a,b))$root

    root<-function(theta)
    {
      sfabz(theta,alpha=alpha,psi=psi) - (1- pt( (theta-y)/s,dof )/alpha)
    }
    a<-b<-y+s*qt(alpha,dof)
    #while(root(a)>0){ a<- a + s*qnorm(alpha)*.25 }
    #while(root(b)<0){ b<- b + s*qnorm(1-alpha)*.25 }
    while(root(a)>0){ a<- a - s }
    while(root(b)<0){ b<- b + s }
    thetaL<-uniroot(root,c(a,b))$root
  }

c(thetaL,thetaU)
}


#' @title Bayes-optimal spending function
#'
#' @description Compute Bayes optimal spending function
#'
#' @details This function computes the 
#' value of s that minimizes the acceptance probability of a 
#' biased level-alpha test for a normal population with 
#' known variance, under a specified  prior
#' predictive distribution.
#' 
#' @param theta value of theta being tested
#' @param psi a list of parameters for the spending function, including 
#' \enumerate{
#' \item mu, the prior expectation of E[y]
#' \item tau2, the prior variance of E[y]
#' \item sigma2 the variance of y
#' }
#' @param alpha level of test
#'
#' @author Peter Hoff 
#' 
#' @export
sfabz<-function(theta,psi,alpha=.05)
{
  mu<-psi$mu ; tau2<-psi$tau2 ; sigma2<-psi$sigma2

  s<-1*(theta>mu)
  if(tau2>0)
  {
    igfun<-function(x,alpha)
    {
    gsmx <-function(s){ qnorm(alpha*s) - qnorm(alpha*(1-s)) - x }
    uniroot(gsmx,interval=c(0,1),maxiter=2000,tol=.Machine$double.eps^0.5)$root
    }
    s<-igfun( 2*sqrt(sigma2)*(theta-mu)/tau2,alpha)
  }
  s
}


