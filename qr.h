// Householder QR algorithm

#include <armadillo>
#include <cmath>

/// @brief Calculates eigenvalues/eigenvecs of a matrix using QR algorithm with Householder reflections
/// @param A  symmetric matrix
/// @param eigval  vector to store eigenvalues
/// @param eigvec  matrix to store eigenvectors
/// @param maxIterations  maximum number of iterations
/// @param tol  tolerance for convergence
/// @param verbose  print convergence information
/// @param method  QR decomposition method: "h" for Householder, "gs" for Gram-Schmidt, "g" for Givens
void qrAlgorithm(arma::mat& A, 
                arma::vec &eigval, 
                arma::mat &eigvec, 
                int maxIterations = 10000, 
                double tol = 1e-6, 
                bool verbose=false, 
                std::string method="h");