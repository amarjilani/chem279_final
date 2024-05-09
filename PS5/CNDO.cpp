#include "CNDO.h"
#include "../qr.h"

// function to calculate [0][o] 
double boys_function(AO &sh1, AO &sh2, int k, int kp, int l, int lp) {

  // calculate all the terms 
  double sigma_a = 1.0 / (sh1.get_alpha()(k) + sh1.get_alpha()(kp));
  double sigma_b = 1.0 / (sh2.get_alpha()(l) + sh2.get_alpha()(lp));
  double Ua = pow((M_PI * sigma_a), 1.5); 
  double Ub = pow((M_PI * sigma_b), 1.5); 
  double Vsq = 1.0 / (sigma_a + sigma_b);
  arma::vec RA = sh1.get_R0(); 
  arma::vec RB = sh2.get_R0(); 
  arma::vec RAB = RA - RB;     
  double Rab_squared = dot(RAB, RAB);
  double T = Vsq * Rab_squared;
  double res; 

  // calculate 
  if (T == 0){
    res =(Ua * Ub) * sqrt(2 * Vsq) * sqrt(2/M_PI);
  } else {
    res = (Ua * Ub) * (1/sqrt(Rab_squared)) * erf(sqrt(T));
  }
  return res; 
}


// calculate the gamma value between two AOs
double gamma(AO &sh1, AO &sh2){
  double g = 0; 
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          double dk = sh1.get_d_coe()(i);
          double dkprime = sh1.get_d_coe()(j);
          double dl = sh2.get_d_coe()(k);
          double dlprime = sh2.get_d_coe()(l);
          double res = dk * dkprime * dl * dlprime * boys_function(sh1, sh2, i, j, k, l);
          g += res;
        }
      }
    }
  }
  // convert to ev 
  return g * 27.2114; 
  
}

// generate the gamma matrix
void generate_GammaMat(std::vector<AO> &MoleculeAOs, arma::mat &Gamma_mat) {
  int dim = 2; 
  for (size_t k = 0; k <dim; k++){
    for (size_t j = 0; j < dim; j++){
      double G_1AO = gamma(MoleculeAOs[k], MoleculeAOs[j]);
      Gamma_mat(k,j) = G_1AO;
    }
  }
}

// calculate the total electron density 
arma::vec P_AA(std::vector<int> atoms, arma::mat &P_alpha, arma::mat &P_beta) {
    
    arma::mat P_total = P_alpha + P_beta;
    arma::vec P_AA = arma::zeros<arma::vec>(atoms.size());
    for (size_t i = 0; i < atoms.size(); i++) {
      if (atoms[i] == 1) {
        P_AA(i) = P_total(i, i);
      } else {
        P_AA(i) = P_total(i, i) + P_total(i + 1, i + 1) + P_total(i+2, i+2) + P_total(i+3, i+3);
      }
    }
    return P_AA;

}

// build the Fock matrix
void build_Fock_matrix(std::vector<AO> &MoleculeAOs, std::vector<int> atoms, arma::mat &S, arma::mat &G_mat, arma::mat &P_alpha, arma::mat &P_beta, arma::mat &F_alpha, int valence_elec) {
    int dim = MoleculeAOs.size();
    F_alpha.zeros(dim, dim);
    arma::mat P_total = P_alpha + P_beta; 
    arma::vec Paa = P_AA(atoms, P_alpha, P_beta);

    for (size_t mu = 0; mu < dim; mu++) {
      // get curent atom index
        int atomIdxA = MoleculeAOs[mu].get_atom_idx(); 
        double Z_A = valence_map.at(MoleculeAOs[mu].get_atomicN());
        double gamma_AA = G_mat(atomIdxA, atomIdxA); 
        double gamma_AB = G_mat(0, 1); 

        for (size_t nu = 0; nu < dim; nu++) {
            if (mu == nu) {
                // extract empiracal data
                double value_mu = (MoleculeAOs[mu].get_lable().find("s") != std::string::npos) ?
                                  elementData.at(MoleculeAOs[mu].get_atomicN()).Is_As :
                                  elementData.at(MoleculeAOs[mu].get_atomicN()).Ip_Ap;

                double sum_pBB_ZB_gamma_AB = 0.0; 

                // sum over all atoms
                for (size_t i = 0; i < atoms.size(); i++) {
                    if (i != atomIdxA) { 
                        sum_pBB_ZB_gamma_AB += (Paa[i] - valence_map.at(atoms[i])) * gamma_AB;
                    }
                }
                // diagonal calculation
                F_alpha(mu, mu) = -value_mu + (Paa[atomIdxA] - Z_A - (P_alpha(mu, mu) - 0.5)) * gamma_AA + sum_pBB_ZB_gamma_AB;
            } else {
                // Off-diagonal elements of F_alpha
                double g_off_diag = G_mat(atomIdxA, MoleculeAOs[nu].get_atom_idx());
                double beta_A = elementData.at(MoleculeAOs[mu].get_atomicN()).beta;
                double beta_B = elementData.at(MoleculeAOs[nu].get_atomicN()).beta;

                F_alpha(mu, nu) = 0.5 * (-beta_A - beta_B) * S(mu, nu) - P_alpha(mu, nu) * g_off_diag;
                F_alpha(nu, mu) = F_alpha(mu, nu);
            }
        }
    }
}

// calculate the nuclear repulsion energy
double calculate_nuclear_repulsion(std::vector<AO> MoleculeAOs, std::vector<int> atoms) {
    double nuclear_repulsion = 0.0;
    // find MO that belongs to atom index 0 to extract pos 
    arma::vec A1_pos;
    for (size_t i = 0; i < MoleculeAOs.size(); i++) {
        if (MoleculeAOs[i].get_atom_idx() == 0) {
            A1_pos = MoleculeAOs[i].get_R0();
            break;
        }
    }
    
    // find MO that belongs to atom index 1 to extract pos.
    arma::vec A2_pos;
    for (size_t i = 0; i < MoleculeAOs.size(); i++) {
        if (MoleculeAOs[i].get_atom_idx() == 1) {
            A2_pos = MoleculeAOs[i].get_R0();
            break;
        }
    }

    // calculate repulsion 
    for (size_t i = 0; i < atoms.size(); i++) {
        for (size_t j = i + 1; j < atoms.size(); j++) {
            double Z1 = valence_map.at(atoms[i]);
            double Z2 = valence_map.at(atoms[j]);
            double distance = sqrt(pow(A1_pos(0) - A2_pos(0), 2) + pow(A1_pos(1) - A2_pos(1), 2) + pow(A1_pos(2) - A2_pos(2), 2));
            nuclear_repulsion += Z1 * Z2 / distance;
        }
    }

    // conver to ev 
    return nuclear_repulsion * 27.2114;
}

// calculate electronic energy using h_core, F_alpha, F_beta, P_alpha, P_beta
double calculate_electronic_energy(const arma::mat& H_core, const arma::mat& F_alpha, const arma::mat& F_beta, const arma::mat& P_alpha, const arma::mat& P_beta) {
    double electronicEnergy = 0.0;
    for (arma::uword mu = 0; mu < P_alpha.n_rows; ++mu) {
        for (arma::uword nu = 0; nu < P_alpha.n_cols; ++nu) {
            electronicEnergy += P_alpha(mu, nu) * (H_core(mu, nu) + F_alpha(mu, nu));
            electronicEnergy += P_beta(mu, nu) * (H_core(mu, nu) + F_beta(mu, nu));
        }
    }
    return electronicEnergy *= 0.5; 

}

// perform SCF to converge the density matrix and calculate energy
void SCF(std::vector<AO> &MoleculeAOs, std::vector<int> atoms, arma::mat &OV_mat, arma::mat &G_mat, arma::mat &H_core, arma::mat &P_alpha, arma::mat &P_beta, arma::mat &F_alpha, arma::mat &F_beta, int p, int q, int valence_elec, std::string method) {
    int dim = MoleculeAOs.size();
    int max_iter = 100;
    double convergence = 1e-6;
    for (size_t iter = 0; iter < max_iter; iter++) {
        // build Fock matrix

        build_Fock_matrix(MoleculeAOs, atoms, OV_mat, G_mat, P_alpha, P_beta, F_alpha, valence_elec);
        build_Fock_matrix(MoleculeAOs, atoms, OV_mat, G_mat, P_beta, P_alpha, F_beta, valence_elec);

        // assumng initial guess of 0 
        if (iter == 0) {
          H_core = F_alpha;
        }

        // solve the eigenvalue problem
        arma::vec eigval_alpha, eigval_beta;
        arma::mat eigvec_alpha, eigvec_beta;

        if (method != "e") {
          qrAlgorithm(F_alpha, eigval_alpha, eigvec_alpha, 10000, 1e-6, false, method);
          qrAlgorithm(F_beta, eigval_beta, eigvec_beta, 10000, 1e-6, false, method);
        } else {
          arma::eig_sym(eigval_alpha, eigvec_alpha, F_alpha);
          arma::eig_sym(eigval_beta, eigvec_beta, F_beta);
        }

        // sort eigen values
        arma::uvec idx_alpha = arma::sort_index(eigval_alpha);
        arma::uvec idx_beta = arma::sort_index(eigval_beta);

        // sort eigenvectors
        eigval_alpha = eigval_alpha(idx_alpha);
        eigvec_alpha = eigvec_alpha.cols(idx_alpha);
        eigval_beta = eigval_beta(idx_beta);
        eigvec_beta = eigvec_beta.cols(idx_beta);


        // print

        // construct new density matrices for alpha and beta electrons
        arma::mat Pa_new = arma::zeros(dim, dim);
        arma::mat Pb_new = arma::zeros(dim, dim);
        
        // populate with only occupied orbitals
        for (int mu = 0; mu < p; mu++) {
            Pa_new += eigvec_alpha.col(mu) * eigvec_alpha.col(mu).t();
        }
        for (int mu = 0; mu < q; mu++) {
            Pb_new += eigvec_beta.col(mu) * eigvec_beta.col(mu).t();
        }

        // print

        // check if Pa_new and Pb_new are converged
        if (arma::approx_equal(P_alpha, Pa_new, "absdiff", convergence) && arma::approx_equal(P_beta, Pb_new, "absdiff", convergence)) {
            // calculate energy 
            double E = 0.0;
            double nuclear_repulsion = calculate_nuclear_repulsion(MoleculeAOs, atoms);
            std::cout << "Nuclear repulsion: " << nuclear_repulsion << " eV" << std::endl;
            double electronicEnergy = calculate_electronic_energy(H_core, F_alpha, F_beta, P_alpha, P_beta);
            std::cout << "Electronic energy: " << electronicEnergy << " eV" << std::endl; 
            std::cout << "Total molecular energy: " << nuclear_repulsion + electronicEnergy << " eV" << std::endl; 

            break;
        } else {
            P_alpha = Pa_new;
            P_beta = Pb_new;
        }

    }
}

/*
------------------------------------------------------------------------------------------------------------------------
PROBLEM SET 5 STARTS HERE
------------------------------------------------------------------------------------------------------------------------
*/ 

// calculate overlap derivative

/// @brief Function to calculate the derivative of the 1D overlap interval
/// @param xa coordinate of the first atom
/// @param xb coordinate of the second atom
/// @param alphaa alpha of the first atom
/// @param alphab alpha of the second atom
/// @param la lmn of the first atom
/// @param lb lmn of the second atom
/// @return 
double overlapDerivative1D(double xa, double xb, double alphaa, double alphab, int la, int lb){
  double first_term = -la * Overlap_onedim(xa, xb, alphaa, alphab, la - 1, lb);
  double second_term = 2 * alphaa * Overlap_onedim(xa, xb, alphaa, alphab, la + 1, lb);
  return first_term + second_term;

}

/// @brief Calculate the derivative of the 3D overlap interval
/// @param sh1 shell 1
/// @param sh2 shell 2
/// @param k 
/// @param l 
/// @return 
arma::vec overlapDerivative3D(AO &sh1, AO &sh2, int k, int l) {
  // 1D overlap
  double Ix = Overlap_onedim(sh1.get_R0()(0), sh2.get_R0()(0), sh1.get_alpha()(k), sh2.get_alpha()(l), sh1.get_lmn()(0), sh2.get_lmn()(0));
  double Iy = Overlap_onedim(sh1.get_R0()(1), sh2.get_R0()(1), sh1.get_alpha()(k), sh2.get_alpha()(l), sh1.get_lmn()(1), sh2.get_lmn()(1));
  double Iz = Overlap_onedim(sh1.get_R0()(2), sh2.get_R0()(2), sh1.get_alpha()(k), sh2.get_alpha()(l), sh1.get_lmn()(2), sh2.get_lmn()(2));

  // overlap derivatives
  double dIx = overlapDerivative1D(sh1.get_R0()(0), sh2.get_R0()(0), sh1.get_alpha()(k), sh2.get_alpha()(l), sh1.get_lmn()(0), sh2.get_lmn()(0));
  double dIy = overlapDerivative1D(sh1.get_R0()(1), sh2.get_R0()(1), sh1.get_alpha()(k), sh2.get_alpha()(l), sh1.get_lmn()(1), sh2.get_lmn()(1));
  double dIz = overlapDerivative1D(sh1.get_R0()(2), sh2.get_R0()(2), sh1.get_alpha()(k), sh2.get_alpha()(l), sh1.get_lmn()(2), sh2.get_lmn()(2));

  arma::vec overlapDeriv3D = {dIx * Iy * Iz, Ix * dIy * Iz, Ix * Iy * dIz};
  return overlapDeriv3D;
}

/// @brief Calculate the contracted overlap derivative
/// @param sh1 shell 1
/// @param sh2 shell 2
/// @return 
arma::vec contractedOverlapDerivative(AO &sh1, AO &sh2) {
  arma::vec overlapDeriv3D = arma::zeros<arma::vec>(3);
  for (int k = 0; k < 3; k++) {
    for (int l = 0; l < 3; l++) {
      arma::vec overlap = overlapDerivative3D(sh1, sh2, k, l);
      overlapDeriv3D += sh1.get_d_coe()(k) * sh2.get_d_coe()(l) * overlap;
    }
  }
  //PRINT
  return overlapDeriv3D;
}

/// @brief Generate the overlap derivative matrix (Suv_RA)
/// @param MoleculeAOs 
/// @param OV_mat 
void overlapDerivativeMatrix(std::vector<AO> &MoleculeAOs, arma::mat &OV_mat) {
  int dim = MoleculeAOs.size();
  for (size_t k = 0; k < dim; k++) {
    for (size_t j = 0; j < dim; j++) {
      if (MoleculeAOs[k].get_atom_idx() == MoleculeAOs[j].get_atom_idx()) {
        OV_mat.col(k * dim + j) = arma::zeros<arma::vec>(3);
        continue;
      } else {
        arma::vec overlapDeriv3D = contractedOverlapDerivative(MoleculeAOs[k], MoleculeAOs[j]);
        OV_mat.col(k * dim + j) += overlapDeriv3D;
      }
    }
  }
}

/// @brief Calculate derivative of the boysfunction 
/// @param sh1 shell1
/// @param sh2 shell 2
/// @param k 
/// @param kp 
/// @param l 
/// @param lp 
/// @return 
arma::vec boysDerivative(AO &sh1, AO &sh2, int k, int kp, int l, int lp) {

  // calculate all the terms 
  double sigma_a = 1.0 / (sh1.get_alpha()(k) + sh1.get_alpha()(kp));
  double sigma_b = 1.0 / (sh2.get_alpha()(l) + sh2.get_alpha()(lp));
  double Ua = pow((M_PI * sigma_a), 1.5); 
  double Ub = pow((M_PI * sigma_b), 1.5); 
  double Vsq = 1.0 / (sigma_a + sigma_b);
  arma::vec RA = sh1.get_R0(); 
  arma::vec RB = sh2.get_R0(); 
  arma::vec RAB = RA - RB;     
  double Rab_squared = dot(RAB, RAB); 
  double Rab_norm = sqrt(Rab_squared); // | Ra - Rb | 
  double T = Vsq * Rab_squared;
  arma::vec res; 

  arma::vec first_term = ((Ua * Ub) * (RAB)) / pow(Rab_norm, 2);
  double second_term = -(erf(sqrt(T))/ Rab_norm) + ((2 * sqrt(Vsq)) / sqrt(M_PI)) * exp(-T);
  res = first_term * second_term;
  
  // 0 if distance is 0 
  if (T == 0 ) {
    res = arma::zeros<arma::vec>(3);
  }

  return res; 
}

/// @brief Calculate the gamma derivative between two AOs
/// @param sh1 shell 1
/// @param sh2 shell 2
/// @return 
arma::vec gammaDerivative(AO &sh1, AO &sh2){
  arma::vec g = arma::zeros<arma::vec>(3); 
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          double dk = sh1.get_d_coe()(i);
          double dkprime = sh1.get_d_coe()(j);
          double dl = sh2.get_d_coe()(k);
          double dlprime = sh2.get_d_coe()(l);
          g += dk * dkprime * dl * dlprime * boysDerivative(sh1, sh2, i, j, k, l);
        }
      }
    }
  }
  // convert to ev 
  return g * 27.2114; 
  
}

/// @brief Generate gamma matrix derivative
/// @param MoleculeAOs vector of AOs
/// @param Gamma_mat gamma derivative matrix to fill 
void generate_GammaDerivativeMat(std::vector<AO> &MoleculeAOs, arma::mat &Gamma_mat) {
  int dim = 2; 
  for (size_t k = 0; k < dim; k++){
    for (size_t j = 0; j < dim; j++){
      arma::vec G_1AO = gammaDerivative(MoleculeAOs[k], MoleculeAOs[j]);
      Gamma_mat.col(k * dim + j) += G_1AO;
    }
  }
}

/// @brief Calculate the nuclear part of the gradient 
/// @param MoleculeAOs vector of AOs
/// @param info info regarding the molecule
/// @param grad mat to fill 
void nuclearGradient(std::vector<AO> &MoleculeAOs, MolInfo info, arma::mat &grad) {
  for (size_t i = 0; i < info.num_Atoms; i++) {
    int atom1 = info.atoms[i];
    int Za = valence_map.at(atom1);
    for (size_t j = 0; j < info.num_Atoms; j++) {
      int atom2 = info.atoms[j];
      if (i == j) {
        continue;
      }
      int Zb = valence_map.at(atom2);

      arma::vec Ra = arma::zeros<arma::vec>(3);
      arma::vec Rb = arma::zeros<arma::vec>(3);

      for (size_t k = 0; k < MoleculeAOs.size(); k++) {
        if (MoleculeAOs[k].get_atom_idx() == i) {
          Ra = MoleculeAOs[k].get_R0();
        } else if (MoleculeAOs[k].get_atom_idx() == j) {
          Rb = MoleculeAOs[k].get_R0();
        }
      }

      arma::vec RAB = Ra - Rb;
      double RAB_norm = sqrt(dot(RAB, RAB));
      double RAB_norm3 = pow(RAB_norm, 3);

      grad.col(i) += -(Za * Zb) * (RAB/RAB_norm3);
      // conver to ev
      grad.col(i) *= 27.2114;
    }
  }
}


/// @brief Generate the X matrix
/// @param MoleculeAOs vector of AOs 
/// @param Pa density matrix alpha
/// @param Pb density matrix beta
/// @param X_mat mat to fill 
void generateXmat(std::vector<AO> &MoleculeAOs, arma::mat Pa, arma::mat Pb, arma::mat &X_mat) {
  arma::mat Ptot = Pa + Pb;
  for (size_t i = 0; i < MoleculeAOs.size(); i++) {
    for (size_t j = 0; j < MoleculeAOs.size(); j++) {
      double betaA = elementData.at(MoleculeAOs[i].get_atomicN()).beta;
      double betaB = elementData.at(MoleculeAOs[j].get_atomicN()).beta;
      X_mat(i, j) = (-betaA - betaB) * Ptot(i, j);
    }
  }
}

/// @brief Generate the Y matrix
/// @param MoleculeAOs vector of AOs
/// @param info mol info
/// @param Pa density matrix alpha
/// @param Pb density matrix beta
/// @param Y_mat mat to fill 
void generateYmat(std::vector<AO> &MoleculeAOs, MolInfo info, arma::mat Pa, arma::mat Pb, arma::mat &Y_mat) {
  arma::vec Pt = P_AA(info.atoms, Pa, Pb);
  for (int i = 0; i < info.num_Atoms; i++) {
    for (int j = 0; j < info.num_Atoms; j++) {
      // terms 
      double Paa = Pt(i);
      double Pbb = Pt(j);
      double ZA = valence_map.at(info.atoms[i]);
      double ZB = valence_map.at(info.atoms[j]);
      double sum = 0;

      // summation 
      for (int k = 0; k < MoleculeAOs.size(); k++) {
        for (int l = 0; l < MoleculeAOs.size(); l++) {
          if (MoleculeAOs[k].get_atom_idx() == i && MoleculeAOs[l].get_atom_idx() == j) {
            sum += Pa(k,l) * Pa(k,l) + Pb(k,l) * Pb(k,l);
          }
        }
      }
      // combine 
      Y_mat(i, j) = Paa*Pbb - ZB*Paa - ZA*Pbb - sum;
    }
  }
}

/// @brief Calculate the electron portion of the gradient 
/// @param MoleculeAOs vector of AOs
/// @param info mol info
/// @param S_deriv Su_RA
/// @param G_deriv Gamma RA
/// @param X_mat x matrix
/// @param Y_mat y matrix
/// @param grad mat to fill 
void electronGradient(std::vector<AO> &MoleculeAOs, MolInfo info, arma::mat &S_deriv, arma::mat &G_deriv, arma::mat &X_mat, arma::mat &Y_mat, arma::mat &grad) {
  arma::mat grad1 = arma::zeros<arma::mat>(3, info.num_Atoms);
  arma::mat grad2 = arma::zeros<arma::mat>(3, info.num_Atoms);
  // first term 
  for (int i = 0; i < MoleculeAOs.size(); i++) {
    for (int j = 0; j < MoleculeAOs.size(); j++) {
      if (i == j) {
        continue;
      } else {
        grad1.col(MoleculeAOs[i].get_atom_idx()) += X_mat(i, j) * S_deriv.col(i * MoleculeAOs.size() + j);
      }
    }
  }
  // second term 
  for (int i = 0; i < info.num_Atoms; i++) {
    for (int j = 0; j < info.num_Atoms; j++) {
      if (i == j) {
        continue;
      } else {
        grad2.col(i) += Y_mat(i, j) * G_deriv.col(i * info.num_Atoms + j);
      }
    }
  }
  // combine 
  grad = grad1 + grad2;
  }