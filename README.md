# frechetAEE

This repository provides an R software for an asymptotical efficient estimator for the Frechet distribution, $FR(\theta_1, \theta_2)$, 
where the PDF of the distribution is as follows:

$$
f(x) = \theta_1 \theta_2 x^{-(\theta_1 + 1)} \exp(-\theta_2 x^{-\theta_1}), \quad x > 0.
$$

Our function provides estimators for $\theta_1$ and $\theta_2$. Our function is implemented using \code{Rcpp} and \code{RcppArmadillo}, to ensure the practicality and the speed. 
We provide a simple example for utilizing the function:
```
library(Rcpp)
library(VGAM)

sourceCpp("frechet_est.cpp")

n = 100
X = rfrechet(n, location = 0)

est = ce_frechet_cpp(X)
```
Then \code{est} provides the estimates for $\theta_1$ and $\theta_2$.