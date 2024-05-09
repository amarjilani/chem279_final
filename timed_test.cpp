#include "qr.h"
#include "utils/timer.h"


int main() {
    // set rand seed
    arma::arma_rng::set_seed(time(NULL));

    arma::mat A(10, 10, arma::fill::randu);

    // make symmetric
    A = A + A.t();


    // eigen of original
    tools::timer t("arma::eig_sym");
    double elapsed = 0;
    for (int i = 0; i < 100; i++) {
        arma::vec eigval;
        arma::mat eigvec;
        t.start();
        arma::eig_sym(eigval, eigvec, A);
        t.stop();
        elapsed += t.elapsed();
    }

    // average time
    std::cout << "Average time for arma::eig_sym: " << elapsed / 100 << std::endl;

    tools::timer t2("Householder");
    double elapsed_h = 0;
    for (int i = 0; i < 100; i++) {
        arma::vec eigval;
        arma::mat eigvec;
        arma::mat A2 = A;
        t2.start();
        qrAlgorithm(A2, eigval, eigvec, 10000, 1e-6, false, "h");
        t2.stop();
        elapsed_h += t2.elapsed();
    }

    // average time
    std::cout << "Average time for householder QR: " << elapsed_h / 100 << std::endl; 

    // given 
    tools::timer t3("Given");
    double elapsed_g = 0;
    for (int i = 0; i < 100; i++) {
        arma::vec eigval;
        arma::mat eigvec;
        arma::mat A2 = A;
        t3.start();
        qrAlgorithm(A2, eigval, eigvec, 10000, 1e-6, false, "g");
        t3.stop();
        elapsed_g += t3.elapsed();
    }

    // average time
    std::cout << "Average time for Given rotation QR: " << elapsed_g / 100 << std::endl; 

    tools::timer t4("Gram-Shmidt");
    double elapsed_gs = 0;
    for (int i = 0; i < 100; i++) {
        arma::vec eigval;
        arma::mat eigvec;
        arma::mat A2 = A;
        t4.start();
        qrAlgorithm(A2, eigval, eigvec, 10000, 1e-6, false, "gs");
        t4.stop();
        elapsed_gs += t4.elapsed();
    }

    // average time
    std::cout << "Average time for Gram-Schmidt QR: " << elapsed_gs / 100 << std::endl; 


    return 0;
}