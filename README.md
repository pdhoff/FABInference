## FABInference

#### Construction of FAB p-values and confidence intervals


[Package website](https://pdhoff.github.io/FABInference/)


#### About
In multiparameter problems, information sharing across parameters 
can be used to improve the power of statistical hypothesis tests, thereby
providing smaller $p$-values and narrower confidence intervals, on average across parameters. 
The `FABInference` package provides information sharing in linear and 
generalized linear regression models using a syntax similar to the 
built-in R functions `lm` and `glm`. 

#### Usage
Suppose you want to get FAB $p$-values for the predictors $x_{i,1},\ldots, x_{i,p}$ in the linear model

\[
 y_i = \alpha_0 + \alpha_1 w_{i,1} + \alpha_2 w_{i,2} + \beta_1 x_{i,1} + \cdots + \beta_p x_{i,p} + \epsilon_i,  
\]

where $w_{i,1}$ and $w_{i,2}$ (and potentially other $w_{i,j}$'s) are additional control variables you'd like to have in the model. Then you need to

1. column-bind the $x$-variables into an $n\times p$ matrix `X`, 
   e.g. `X<-cbind(x1,x2,x3)`; 

2. run the command `fit<-lmFAB(y~w1+w2,X)`.  

The output is similar to the output of the `lm` command, so you can 
type `summary(fit)` to see the FAB $p$-values. The FAB $p$-values and 
confidence intervals are stored in `fit$FABpv` and `fit$FABci`. 

If $\beta_1,\ldots, \beta_p$ correspond to $p$ objects about which you 
have additional covariate information (say attributes $\{(v_{j,1},v_{j,2}), 
j =1,\ldots, p\}$ you might be interested in fitting the model 
`fit<-lmFAB(y~w1+w2,X,~v1+v2)`, where `v1` and `v2` are $p$-dimensional 
vectors giving the attributes associated with $\beta_1,\ldots, \beta_p$. 
The additional term specifies a *linking model* for $\beta_1,\ldots, \beta_p$. 
Importantly, the linking model doesn't have to be correct in any way for the 
FAB $p$-values of confidence intervals to be valid. However, the better the linking model, the smaller the $p$-values and the narrower the intervals. 

FAB inference for generalized linear models can be obtained similarly 
using the command `glmFAB`. In this case, the $p$-values and confidence 
intervals are valid asymptotically (just like the standard 
$p$-values and intervals). Fitting a normal linear regression 
with `glmFAB` is much faster than using `lmFAB` because the former uses 
an asymptotic approximation. 

#### Theoretical details
In the simplest case of a normally distributed estimator $\hat\theta$ of
$\theta$ such that $\hat \theta \sim N(\theta,\sigma^2)$, a standard $p$-value
and confidence interval are based on the test statistic $|\hat\theta|$. 
A FAB $p$-value and confidence interval is  based on the statistic
$|\hat\theta + a|$, where $a$ is determined from indirect information about 
the sign and magnitude of $\theta$. The functional form of the FAB $p$-value
is extremely simple: 

\[
 p_{FAB}(\hat\theta,a) = 1- | \Phi(\hat\theta+2a) - \Phi(-\hat\theta) |, 
\]

where $\Phi$ is the standard normal CDF. The FAB confidence interval is a bit more complicated. In multiparameter settings, the optimal choice for $a$ for one parameter may be estimated from data on the other parameters, using a linking model that relates the parameters to each other. Importantly, the FAB confidence intervals and $p$-values have correct frequentist error rates, even if the linking model is incorrect. 

#### Installation


```{r,eval=FALSE}
# Release version on CRAN
install.packages("FABInference")


# Development version on GitHub 
devtools::install_github("pdhoff/FABInference")  
```


#### Documentation and citation


"Smaller p-values via indirect information". P.D. Hoff.  [arXiv:1907.12589](https://arxiv.org/abs/1907.12589) 

"Exact adaptive confidence intervals for small areas". K. Burris and P.D. Hoff. 
[arXiv:1809.09159](https://arxiv.org/abs/1809.09159) Journal of Survey Statistics and Methodology.

"Exact adaptive confidence intervals for linear regression coefficients". 
P.D. Hoff and C. Yu. [arXiv:1705.08331](https://arxiv.org/abs/1705.08331) Electron. J. Stat., 13(1):94–119, 2019. 

"Adaptive multigroup confidence intervals with constant coverage". 
[arXiv:1612.08287](https://arxiv.org/abs/1612.08287) C. Yu and P.D. Hoff. Biometrika, 105(2):319–335, 2018. 


#### Some examples

* [Demo and simulation study](articles/simstudy.html)

* [Small area estimation](articles/exampleFHmodel.html) Replication file for Hoff(2019)

* [Hidden Markov model](articles/exampleHMM.html) Replication file for Hoff(2019)

* [Linear model interactions](articles/exampleInteraction.html) Replication file for Hoff(2019)

* [Logistic regression](articles/exampleLogistic.html) Replication file for Hoff(2019)


