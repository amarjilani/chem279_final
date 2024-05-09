#include <iostream>
#include <armadillo>
#include <chrono>
#include "qr.h"


/// @brief Generates a symmetric sparse matrix
/// @param size Matrix size
/// @param sparsity Fraction of non-zero elements (0.0 to 1.0)
arma::mat generateSparseSymmetricMatrix(int size, double sparsity) {
    arma::sp_mat sparseMat = arma::sprandu<arma::sp_mat>(size, size, sparsity);
    arma::mat denseMat = arma::mat(sparseMat);
    return denseMat + denseMat.t();  // Make symmetric
}

/// @brief Measures the execution time of a function
/// @param func Function to measure
/// @return Duration in seconds
template<typename Func>
double measureExecutionTime(Func&& func) {
    auto start = std::chrono::high_resolution_clock::now();
    func();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    return duration.count();
}

int main(int argc, char** argv) {
    int matrixSize = 10;
    double sparsity = 0.1;

    if (argc >= 2) {
        matrixSize = std::stoi(argv[1]);
    }
    if (argc >= 3) {
        sparsity = std::stod(argv[2]);
    }

    std::cout << "Matrix Size: " << matrixSize << ", Sparsity: " << sparsity << std::endl;

    // generate a sparse symmetric matrix
    arma::mat sparseMat = generateSparseSymmetricMatrix(matrixSize, sparsity);

    arma::vec eigval_h, eigval_g;
    arma::mat eigvec_h, eigvec_g;

    // measure execution time for Householder method
    double time_h = measureExecutionTime([&]() {
        qrAlgorithm(sparseMat, eigval_h, eigvec_h, 1000, 1e-6, false, "h");
    });
    std::cout << "Householder Time: " << time_h << " seconds" << std::endl;

    // measure execution time for Givens method
    double time_g = measureExecutionTime([&]() {
        qrAlgorithm(sparseMat, eigval_g, eigvec_g, 1000, 1e-6, false, "g");
    });
    std::cout << "Givens Time: " << time_g << " seconds" << std::endl;


    return 0;
}
