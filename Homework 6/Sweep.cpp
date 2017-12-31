/*
####################################################
## Stat 202A - Homework 6
## Author: Anoosha Sagar
## Date :
## Description: This script implements QR and Sweep
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


/* ~~~~~~~~~~~~~~~~~~~~~~~~~
   Problem 2: Sweep operator
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

