#include "qr.h"

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

/// @brief Performs Householder QR decomposition on a matrix A
/// @param A matrix to decompose
/// @param Q orthogonal matrix
/// @param R upper triangular matrix
void householderQR(arma::mat& A, arma::mat& Q, arma::mat& R) {
    int n = A.n_rows;
    Q = arma::eye<arma::mat>(n, n); 
    R = A;                           

    for (int k = 0; k < n - 1; k++) {
        arma::vec x = R.submat(k, k, n - 1, k);
        double alpha; 
        if (x(0) >= 0) {
            alpha = -arma::norm(x);
        } else {
            alpha = arma::norm(x);
        }
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

/// @brief Performs Givens rotation on a 2x2 matrix
/// @param a element of matrix
/// @param b element of matrix
/// @param c cosine of rotation angle
/// @param s sine of rotation angle
void givensRotation(double a, double b, double &c, double &s) {
    if (b == 0) {
        c = 1;
        s = 0;
    } else {
        if (abs(b) > abs(a)) {
            double t = -a / b;
            s = 1 / sqrt(1 + t*t);
            c = s*t;
        } else {
            double t = -b / a;
            c = 1 / sqrt(1 + t*t);
            s = c*t;
        }
    }
}

/// @brief Performs QR decomposition using Givens rotations
/// @param A matrix to decompose
/// @param Q orthogonal matrix
/// @param R upper triangular matrix
void givensQR(const arma::mat& A, arma::mat& Q, arma::mat& R) {
    int n = A.n_rows;
    int m = A.n_cols;
    Q = arma::eye<arma::mat>(n, n);
    R = A;

    for (int j = 0; j < m; j++) {
        for (int i = n - 1; i > j; i--) {
            double c, s;
            givensRotation(R(i-1, j), R(i, j), c, s);

            arma::mat G = arma::eye<arma::mat>(n, n);
            G(i-1, i-1) = c;
            G(i, i) = c;
            G(i-1, i) = s;
            G(i, i-1) = -s;

            R = G.t() * R;
            Q = Q * G;
        }
    }
}

/// @brief Performs classical Gram-Schmidt QR decomposition on a matrix A
/// @param A matrix to decompose
/// @param Q orthogonal matrix
/// @param R upper triangular matrix
void gramSchmidtQR(const arma::mat& A, arma::mat& Q, arma::mat& R) {
    int n = A.n_cols;
    Q = arma::mat(A.n_rows, n);
    R = arma::mat(n, n, arma::fill::zeros);

    for (int k = 0; k < n; k++) {
        arma::vec ak = A.col(k);
        for (int j = 0; j < k; j++) {
            R(j, k) = arma::dot(Q.col(j), ak);
            ak -= Q.col(j) * R(j, k);
        }
        R(k, k) = arma::norm(ak);
        Q.col(k) = ak / R(k, k);
    }
}

void qrAlgorithm(arma::mat& A, 
                arma::vec &eigval, 
                arma::mat &eigvec, 
                int maxIterations, 
                double tol, 
                bool verbose, 
                std::string method) {
    int n = A.n_rows;
    arma::mat Q, R;
    arma::mat A0 = A;
    eigvec = arma::eye(n, n);

    for (int iter = 0; iter < maxIterations; iter++) {
        if (method == "h") {
            householderQR(A, Q, R);
        } else if (method == "gs") {
            gramSchmidtQR(A, Q, R);
        } else if (method == "g") {
            givensQR(A, Q, R);
        } else {
            std::cout << "Please select a valid QR decomp method!" << std::endl; 
            return; 
        }

        A = R * Q;
        eigvec *= Q; 

        if (arma::norm(A0 - A, "fro") < tol) { 
            if (verbose) std::cout << "Converged after " << iter << " iterations\n";
            break;
        }
        A0 = A; 
    }
    eigval = A.diag(); 
}




