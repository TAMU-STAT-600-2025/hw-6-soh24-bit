#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Soft-thresholding function, returns scalar
// [[Rcpp::export]]
double soft_c(double a, double lambda){
  // Your function code goes here
  return arma::sign(a) * std::max(std::abs(a) - lambda, 0.0);
}

// Lasso objective function, returns scalar
// [[Rcpp::export]]
double lasso_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& beta, double lambda){
  // Your function code goes here
  
  int n = Xtilde.n_rows;
  
  // Define Residual
  arma::colvec r = Ytilde - Xtilde * beta;
  
  // Same objective function as in R version
  double fobj = (1.0 / (2.0 * n)) * arma::dot(r, r) + lambda * arma::accu(arma::abs(beta));
  
  return fobj;
}

// Lasso coordinate-descent on standardized data with one lamdba. Returns a vector beta.
// [[Rcpp::export]]
arma::colvec fitLASSOstandardized_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, double lambda, const arma::colvec& beta_start, double eps = 0.001){
  // Your function code goes here
  
  int n = Xtilde.n_rows;
  int p = Xtilde.n_cols;
  
  // Initialize beta and current objective value
  arma::colvec beta = beta_start;
  double fcurrent = lasso_c(Xtilde, Ytilde, beta, lambda);
  
  arma::colvec resid = Ytilde - Xtilde * beta;
  
  while (true) {
    
    double fprevious = fcurrent;
    
    for (int i = 0; i < p; i++) {
      
      // Define old beta
      double beta_old_i = beta(i);
      
      // Update a and Beta
      double a = (1.0 / n) * arma::dot(Xtilde.col(i), resid) + beta_old_i;
      beta(i) = soft_c(a, lambda);
      
      // Update resid
      if (beta(i) != beta_old_i) {
        resid -= Xtilde.col(i) * (beta(i) - beta_old_i);
      }
    }
    
    fcurrent = lasso_c(Xtilde, Ytilde, beta, lambda);
    
    // Check if we have reached the optimum
    if (std::abs(fprevious - fcurrent) < eps) {
      break;
    }
    
  }
  
  return beta;
  
}  

// Lasso coordinate-descent on standardized data with supplied lambda_seq. 
// You can assume that the supplied lambda_seq is already sorted from largest to smallest, and has no negative values.
// Returns a matrix beta (p by number of lambdas in the sequence)
// [[Rcpp::export]]
arma::mat fitLASSOstandardized_seq_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& lambda_seq, double eps = 0.001){
  // Your function code goes here
  
  int p = Xtilde.n_cols;
  int n_lambda = lambda_seq.n_elem;
  
  // Initialize a beta matrix
  arma::mat beta_mat(p, n_lambda);
  
  // Warm start
  arma::colvec beta_current(p, arma::fill::zeros);
  
  for (int i = 0; i < n_lambda; i++) {
    
    // Update beta_mat with beta_matrix with corresponding lambda values
    arma::colvec currentfit = fitLASSOstandardized_c(Xtilde, Ytilde, lambda_seq(i), beta_current, eps);
    beta_mat.col(i) = currentfit;
    
    beta_current = currentfit;
  }
  return beta_mat;
}