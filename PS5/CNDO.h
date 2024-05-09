#include "AO.h"

void generate_GammaMat(std::vector<AO> &MoleculeAOs, arma::mat &Gamma_mat);

void build_Fock_matrix(std::vector<AO> &MoleculeAOs, std::vector<int> atoms, arma::mat &S, arma::mat &G_mat, arma::mat &P_alpha, arma::mat &P_beta, arma::mat &F_alpha, int valence_elec);

void SCF(std::vector<AO> &MoleculeAOs, std::vector<int> atoms,  arma::mat &OV_mat, arma::mat &G_mat, arma::mat &H_core, arma::mat &P_alpha, arma::mat &P_beta, arma::mat &F_alpha, arma::mat &F_beta, int p, int q, int valence_elec, std::string method);

void overlapDerivativeMatrix(std::vector<AO> &MoleculeAOs, arma::mat &OV_mat);

void generate_GammaDerivativeMat(std::vector<AO> &MoleculeAOs, arma::mat &Gamma_mat);

void nuclearGradient(std::vector<AO> &MoleculeAOs, MolInfo info, arma::mat &grad);

void generateXmat(std::vector<AO> &MoleculeAOs, arma::mat Pa, arma::mat Pb, arma::mat &X_mat);

void generateYmat(std::vector<AO> &MoleculeAOs, MolInfo info, arma::mat Pa, arma::mat Pb, arma::mat &Y_mat);

void electronGradient(std::vector<AO> &MoleculeAOs, MolInfo info, arma::mat &S_deriv, arma::mat &G_deriv, arma::mat &X_mat, arma::mat &Y_mat, arma::mat &grad);