#' Calculating an aysmptotical efficient estimate for the Frechet distribution
#'
#' @param x A numeric vector.
#' @return
#' \item{theta1}{an estimated shape parameter}
#' \item{theta2}{an estimated scale parameter}
#' \item{l_n_deriv}{a vector containing the first derivative of log likelihood function}
#' \item{fisher_1}{a matrix of the expected fisher information divided by the number of observations}
#' @name ce_frechet_cpp
#' @useDynLib frechetAEE, .registration = TRUE
#' @import Rcpp
#' @export
ce_frechet_cpp
