# frechetAEE

This repository provides an R software for an asymptotical efficient estimator for the Frechet distribution, $FR(\theta_1, \theta_2)$, 
where the PDF of the distribution is as follows:
$$
f(x) = \theta_1 \theta_2 x^{-(\theta_1 + 1)} \exp(-\theta_2 x^{-\theta_1}), \quad x > 0.
$$

Our functions are implemented using \code{Rcpp}, to ensure the practicality and the speed. 

