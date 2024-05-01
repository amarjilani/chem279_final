#include "householder.h"

/// @brief Performs Householder tridiagonalization on a symmetric matrix A inplace
/// @param A  symmetric matrix
void householder(arma::mat& A) {
    int n = A.n_rows;
    
    for (int k = 0; k < n - 2; k++) {        
        arma::vec x = A.submat(k + 1, k, n - 1, k);
        double alpha;   
        if (x(0) >= 0) {
            alpha = -arma::norm(x);
        } else {
            alpha = arma::norm(x);
        }
        arma::vec u = x;
        u(0) -= alpha;
        arma::vec v = u / arma::norm(u);

        arma::mat P = arma::eye<arma::mat>(n - k - 1, n - k - 1) - 2 * v * v.t();
        arma::mat H = arma::eye<arma::mat>(n, n);
        H.submat(k + 1, k + 1, n - 1, n - 1) = P;

        A = H.t() * A * H;
    }
}

void householderQR(arma::mat& A, arma::mat& Q, arma::mat& R) {
    int n = A.n_rows;
    Q = arma::eye<arma::mat>(n, n); 
    R = A;                           

    for (int k = 0; k < n - 1; k++) {
        arma::vec x = R.submat(k, k, n - 1, k);
        double alpha = -arma::norm(x) * (x(0) >= 0 ? 1 : -1);
        arma::vec u = x;
        u(0) -= alpha;
        arma::vec v = u / arma::norm(u);

        arma::mat P = arma::eye<arma::mat>(n - k, n - k) - 2 * v * v.t();
        arma::mat H = arma::eye<arma::mat>(n, n);
        H.submat(k, k, n - 1, n - 1) = P;

        R = H * R; 
        Q = Q * H;
    }
}


void qrAlgorithm(arma::mat& A, int maxIterations = 10000, double tol = 1e-6) {
    int n = A.n_rows;
    arma::mat Q, R;
    arma::mat A0 = A;

    for (int iter = 0; iter < maxIterations; iter++) {
        householderQR(A, Q, R);
        A = R * Q;
        bool converged = true;
        if (arma::norm(A0 - A) < tol) { 
            std::cout << "Converged after " << iter << " iterations\n";
            break;
        }
        A0 = A; 
    }
}

int main() {
    // set rand seed
    arma::arma_rng::set_seed(time(NULL));

    arma::mat A(10, 10, arma::fill::randu);

    // make symmetric
    A = A + A.t();


    std::cout << "Original Matrix A:\n" << A << std::endl;

    // eigen of original
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, A);
    std::cout << "Eigenvalues (A):\n" << eigval << std::endl;

    qrAlgorithm(A);
    std::cout << "After QR Algorithm:\n" << A << std::endl;

    std::cout << "Eigenvalues (Diagonal of A after QR Algorithm):\n" << A.diag() << std::endl;

    return 0;
}

// try out given rotation 

