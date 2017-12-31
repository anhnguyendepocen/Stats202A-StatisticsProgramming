/*
####################################################
## Stat 202A - Homework 2
## Author: 
## Date : 
## Description: This script implements linear regression 
## using the sweep operator
####################################################
 
###########################################################
## INSTRUCTIONS: Please fill in the missing lines of code
## only where specified. Do not change function names, 
## function inputs or outputs. MAKE SURE TO COMMENT OUT ALL 
## OF YOUR EXAMPLES BEFORE SUBMITTING.
##
## Very important: Do not change your working directory
## anywhere inside of your code. If you do, I will be unable 
## to grade your work since R will attempt to change my 
## working directory to one that does not exist.
###########################################################
 
*/ 


# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;



/* ~~~~~~~~~~~~~~~~~~~~~~~~~ 
   Problem 1: Sweep operator 
   ~~~~~~~~~~~~~~~~~~~~~~~~~ */

// [[Rcpp::export()]]
mat mySweepC(const mat A, int m){
  
  /*
  Perform a SWEEP operation on A with the pivot element A[m,m].
  
  A: a square matrix (mat).
  m: the pivot element is A[m, m]. 
  Returns a swept matrix B (which is m by m).
  
  Note the "const" in front of mat A; this is so you
  don't accidentally change A inside your code.
  
  #############################################
  ## FILL IN THE BODY OF THIS FUNCTION BELOW ##
  #############################################
  */
  
  mat B = A;
  int n = B.n_rows;
  int k, i, j;
  
  for(k = 0; k < m; k++){
    
    for( i = 0; i< n; i++){
      for(j = 0; j < n; j++){
        if(i != k && j!= k){
          B(i, j) = B(i, j) - B(i, k) * B(k, j) / B(k, k);
        }
      }
    }
    
    for(i = 0; i < n; i++){
      if(i != k){
        B(i, k) = B(i, k) / B(k, k);
      }
    }
    
    for(j = 0; j < n; j++){
      if(j != k){
        B(k, j) = B(k, j) / B(k, k);
      }
    }
    
    B(k, k) = (-1.0)/B(k, k);
  }
  
  

  
  // Return swept matrix B
  return(B);
    
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
 Problem 2: Linear regression using the sweep operator 
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


// [[Rcpp::export()]]
mat myLinearRegressionC(const mat X, const mat Y){
  
  /*  
  Find the regression coefficient estimates beta_hat
  corresponding to the model Y = X * beta + epsilon
  Your code must use the sweep operator you coded above.
  Note: we do not know what beta is. We are only 
  given a matrix X and a matrix Y and we must come 
  up with an estimate beta_hat.
  
  X: an 'n row' by 'p column' matrix of input variables.
  Y: an 'n row' by '1 column' matrix of responses
    
  #############################################
  ## FILL IN THE BODY OF THIS FUNCTION BELOW ##
  #############################################
  */  
  
  // Let me start things off for you...
  int n = X.n_rows;
  int p = X.n_cols;
  
  mat z(n, 1+p+Y.n_cols);
  mat A, B, beta_hat;
  
  z = join_rows(ones<mat>(1, n).t(), X);
  z = join_rows(z, Y);
  A = z.t() * z;
  
  B = mySweepC(A, p+1);
  beta_hat = B.submat(0, p+1, p, p+1);
  
  
  // Function returns the 'p+1' by '1' matrix 
  // beta_hat of regression coefficient estimates
  return(beta_hat);
      
}
  
  
  
