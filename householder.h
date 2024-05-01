// Householder QR algorithm

#include <armadillo>
#include <cmath>

/// @brief Calculates eigenvalues of a matrix using QR algorithm with Householder reflections
void qrAlgorithm(arma::mat& A, int maxIterations = 10000, double tol = 1e-6);