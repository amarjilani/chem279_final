#include "qr.h"


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

    // print eigen values 
    std::cout << "Eigenvalues (Armadillo):\n" << eigval << std::endl;
    //print eigen vectors 
    std::cout << "Eigenvectors (Armadillo):\n" << eigvec << std::endl;

    arma::vec eigval2;
    arma::mat eigvec2;
    arma::mat A2 = A;
    qrAlgorithm(A2, eigval2, eigvec2, 10000, 1e-6, true, "h");
    
    std::cout << "After QR Algorithm:\n" << A2 << std::endl;
    // print eigen vec
    // sort eigenvecs by eigenvalues
    eigvec2 = eigvec2.cols(arma::sort_index(eigval2));
    std::cout << "Eigenvectors (Householder):\n" << eigvec2 << std::endl;

    eigval2 = arma::sort(eigval2);   
    std::cout << "Eigenvalues (Householder):\n" << eigval2 << std::endl;

    // given
    arma::vec eigval3;
    arma::mat eigvec3;
    arma::mat A3 = A;
    qrAlgorithm(A3, eigval3, eigvec3, 10000, 1e-6, true, "g");

    // print eigen vec
    // sort eigenvecs by eigenvalues
    eigvec3 = eigvec3.cols(arma::sort_index(eigval3));
    std::cout << "Eigenvectors (Givens):\n" << eigvec3 << std::endl;

    eigval3 = arma::sort(eigval3);
    std::cout << "Eigenvalues (Givens):\n" << eigval3 << std::endl;

    // gram schmidt
    arma::vec eigval4;
    arma::mat eigvec4;
    arma::mat A4 = A;
    qrAlgorithm(A4, eigval4, eigvec4, 10000, 1e-6, true, "gs");

    // print eigen vec
    // sort eigenvecs by eigenvalues
    eigvec4 = eigvec4.cols(arma::sort_index(eigval4));
    std::cout << "Eigenvectors (Gram-Schmidt):\n" << eigvec4 << std::endl;

    eigval4 = arma::sort(eigval4);
    std::cout << "Eigenvalues (Gram-Schmidt):\n" << eigval4 << std::endl;


    return 0;
}