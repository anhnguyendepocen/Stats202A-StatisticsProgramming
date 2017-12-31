/*
#########################################################
## Stat 202A - Homework 3
## Author: 
## Date : 
## Description: This script implements QR decomposition,
## linear regression, and eigen decomposition / PCA 
## based on QR.
#########################################################
 
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
 Sign function for later use 
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

// [[Rcpp::export()]]
double signC(double d){
  return d<0?-1:d>0? 1:0;
}



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   Problem 1: QR decomposition 
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */  
  

// [[Rcpp::export()]]
List myQRC(const mat A){ 
  
  /*
  Perform QR decomposition on the matrix A
  Input: 
  A, an n x m matrix (mat)

  #############################################
  ## FILL IN THE BODY OF THIS FUNCTION BELOW ##
  #############################################
  
  */ 
  int row  = A.n_rows;
  int column = A.n_cols;
  mat R = mat(A);
  mat Q = eye<mat>(row, row);
  int k, i;
  
  for(k = 0; k < (column - 1); k++){
    mat x = zeros<mat>(row, 1);
    for(i = k; i < row; i++){
      x(i, 0) = R(i, k);
    }
    mat v = x;
    v(k) = x(k) + signC(x(k, 0)) * norm(x,"fro");
    double s = norm(v, "fro");
    
    if( s != 0){
      mat u = v / s;
      R = R - 2 * (u * (u.t() * R));
      Q = Q - 2 * (u * (u.t() * Q));
    }
    
  }
  
  List output;
  
  
  // Function should output a List 'output', with 
  // Q.transpose and R
  // Q is an orthogonal n x n matrix
  // R is an upper triangular n x m matrix
  // Q and R satisfy the equation: A = Q %*% R
  output["Q"] = Q.t();
  output["R"] = R;
  return(output);
  

}
  
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   Problem 2: Linear regression using QR 
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  
  
// [[Rcpp::export()]]
mat myLinearRegressionC(const mat X, const mat Y){
    
  /*  
  Perform the linear regression of Y on X
  Input: 
  X is an n x p matrix of explanatory variables
  Y is an n dimensional vector of responses
  Do NOT simulate data in this function. n and p
  should be determined by X.
  Use myQRC inside of this function
  
  #############################################
  ## FILL IN THE BODY OF THIS FUNCTION BELOW ##
  #############################################
  
  */  
  int row = X.n_rows;
  int column = X.n_cols;
  
  mat z(row, 1+column+Y.n_cols);
  mat A, B, beta_ls;
  
  z = join_rows(ones<mat>(1, row).t(), X);
  z = join_rows(z, Y);
  
  List sol = myQRC(z);
  mat R = sol["R"];
  
  mat R1 = R.submat(0, 0, column, column);
  mat Y1 = R.submat(0, column+1, column, column+1);
  beta_ls = solve(R1, Y1);
  
  // Function returns the 'p+1' by '1' matrix 
  // beta_ls of regression coefficient estimates
  return(beta_ls.t());
  
}  

/* ~~~~~~~~~~~~~~~~~~~~~~~~ 
 Problem 3: PCA based on QR 
 ~~~~~~~~~~~~~~~~~~~~~~~~~~ */


// [[Rcpp::export()]]

List myEigen_QRC(const mat A, const int numIter = 1000){
  
  /*  
  
  Perform PCA on matrix A using your QR function, myQRC.
  Input:
  A: Square matrix
  numIter: Number of iterations
   
  #############################################
  ## FILL IN THE BODY OF THIS FUNCTION BELOW ##
  #############################################
   
   */  
  int rows = A.n_rows;
  int columns  = A.n_cols;
  mat v = randu<mat>(rows, rows);
  mat temp = v;
  int k;
  List QR;
  
  for(k = 0; k < numIter; k++ ){
  QR = myQRC(v);
  mat q = QR["Q"];
  mat r = QR["R"];
   v = A * q;
   }
  List output_qr = myQRC(v);
  mat Q = output_qr["Q"];
  mat R = output_qr["R"];
  mat rDiag = R.diag();
   List output;
  // // Function should output a list with D and V
  // // D is a vector of eigenvalues of A
  // // V is the matrix of eigenvectors of A (in the
  // // same order as the eigenvalues in D.)
  rDiag.reshape(1,rDiag.n_rows);
  output["D"] = rDiag;
  output["V"] = Q;
  //List output;
  return(output);

}
