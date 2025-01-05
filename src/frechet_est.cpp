#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
Rcpp::List consistent_frechet_cpp(const arma::vec& x) {
  int n = x.n_elem;
  double pi = M_PI;

  arma::vec logx = arma::log(x);
  double squared = arma::sum(arma::square(logx));
  double mean = arma::mean(logx);
  double gamma = -R::digamma(1);

  double theta1_est = pi / std::sqrt( 6* (squared/n - std::pow(mean, 2)));
  double theta2_est = std::exp(theta1_est*mean - gamma);

  return Rcpp::List::create(
    Rcpp::Named("theta1_est") = theta1_est,
    Rcpp::Named("theta2_est") = theta2_est
  );
}

// [[Rcpp::export]]
Rcpp::List ce_frechet_cpp(const arma::vec& x) {
  int n = x.n_elem;
  double pi = M_PI;

  arma::vec logx = arma::log(x);
  double squared = arma::sum(arma::square(logx));
  double mean = arma::mean(logx);
  double gamma = -R::digamma(1);

  double theta1_est_temp = pi / std::sqrt( 6* (squared/n - std::pow(mean, 2)));
  double theta2_est_temp = std::exp(theta1_est_temp*mean - gamma);
  arma::vec xtheta1 = arma::pow(x,theta1_est_temp);

  arma::vec l_n_deriv = {n/theta1_est_temp - arma::sum(logx) + theta2_est_temp*arma::sum(logx / xtheta1),
                         n/theta2_est_temp - arma::sum(1/xtheta1)};
  arma::mat fisher_1(2, 2, arma::fill::zeros);
  fisher_1(0, 0) = (std::pow(theta1_est_temp, 2) * 6) / std::pow(pi, 2);
  fisher_1(0, 1) = fisher_1(1, 0) = (theta1_est_temp * theta2_est_temp * (std::log(theta2_est_temp) + gamma - 1) * 6) / std::pow(pi, 2);
  fisher_1(1, 1) = (std::pow(theta2_est_temp, 2) *
    (std::pow(std::log(theta2_est_temp) + gamma, 2) -
    2 * (std::log(theta2_est_temp) + gamma) +
    1 + std::pow(pi, 2) / 6)* 6) / std::pow(pi, 2);

  arma::vec resid = (fisher_1*l_n_deriv)/n;
  double theta1_est = theta1_est_temp + resid[0];
  double theta2_est = theta2_est_temp + resid[1];
  //arma::mat fisher_1(2,2);
  //fisher_1(0,0) = std::pow(theta1_est_temp, 2);
  //fisher_1(0,1) = theta1_est_temp*theta2_est_temp*(std::log(theta2_est_temp) + gamma - 1);
  //fisher_1(1,0) = theta1_est_temp*theta2_est_temp*(std::log(theta2_est_temp) + gamma - 1);
  //fisher_1(1,1) = std::pow(theta2_est_temp, 2)*( std::pow(std::log(theta2_est_temp) + gamma, 2) - 2*(std::log(theta2_est_temp) + gamma) + 1 + std::pow(pi, 2)/6);

  return Rcpp::List::create(
    Rcpp::Named("theta1_est") = theta1_est,
    Rcpp::Named("theta2_est") = theta2_est,
    Rcpp::Named("l_n_deriv") = l_n_deriv,
    Rcpp::Named("fisher_1") = fisher_1
  );
}


// [[Rcpp::export]]
double sumAbsPairwiseDiff(const arma::vec& x) {
  // Initialize sum
  double result = 0.0;
  arma::uword n = x.n_elem;

  // Pairwise absolute difference sum
  for (arma::uword i = 0; i < n; ++i) {
    if (i % 10000 == 0 || i == n - 1) {
      double percentage = (100.0 * i) / n;
      Rcpp::Rcout << "Iteration: " << i << " / " << n
                  << " (" << percentage << "% completed)" << std::endl;
    }

    for (arma::uword j = 0; j < n; ++j) {
      result += std::abs(x[i] - x[j]);
    }
  }

  // Return the total sum
  return result;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//
