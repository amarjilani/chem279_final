#include <stdlib.h>
#include <stdexcept>
#include <stdio.h>
#include <armadillo>
#include <vector>
#include "AO.h"
#include "CNDO.h"
#include "../utils/timer.h"

using namespace std;

int main(int argc, char* argv[])
{
  arma::mat H_STO3G, C_STO3G, F_STO3G, N_STO3G, O_STO3G;
  H_STO3G.load("basis/H_STO3G.txt");
  C_STO3G.load("basis/C_STO3G.txt");
    F_STO3G.load("basis/F_STO3G.txt");
    N_STO3G.load("basis/N_STO3G.txt");
    O_STO3G.load("basis/O_STO3G.txt");


  vector<AO> MoleculeAOs;
MolInfo info; 
  
  if (argc !=3)
  {
  printf("usage ./hw5 filename eigen_type, for example hw5 example.txt h, eigen_type can be h, g, gs, or e for using armadillo\n");
  return EXIT_FAILURE;
  }
  string fname(argv[1]);
  try
  {
    info = GenerateAOs(MoleculeAOs, fname, H_STO3G, C_STO3G, N_STO3G, O_STO3G, F_STO3G);
  }
  catch (invalid_argument &e)
  {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }

  string method(argv[2]);

  //  gamma 
  arma::mat gamma; 
  gamma.resize(info.atoms.size(), info.atoms.size());
  generate_GammaMat(MoleculeAOs, gamma);


  int dim = MoleculeAOs.size();
  arma::mat OV_mat(dim, dim);
  Eval_OV_mat(MoleculeAOs, OV_mat);


  arma::mat H_core;
  H_core.resize(dim, dim);



  //build Fock matrix
  arma::mat Focka;
  Focka.resize(dim, dim); 
  Focka.zeros();
  arma::mat Fockb;
  Fockb.resize(dim, dim); 
  Fockb.zeros();
  arma::mat Pa(dim, dim);
  arma::mat Pb(dim, dim);
  Pa.zeros();
  Pb.zeros();

  // get energy and density matrix from SCF
  tools::timer t("SCF");
  t.silence();
  double time = 0;
  t.start();
  SCF(MoleculeAOs, info.atoms, OV_mat, gamma, H_core, Pa, Pb, Focka, Fockb, info.p, info.q, info.valence_elec, method);
  t.stop();
  time = t.elapsed();

  // calculate overlap deriv
  arma::mat S_deriv;
  S_deriv.resize(3, dim * dim);
  S_deriv.zeros();
  overlapDerivativeMatrix(MoleculeAOs, S_deriv);
  S_deriv.print("Suv_RA: ");

  // calculate gamma deriv
  arma::mat gamma_deriv;
  gamma_deriv.resize(3, info.atoms.size() * info.atoms.size());
  gamma_deriv.zeros();
  generate_GammaDerivativeMat(MoleculeAOs, gamma_deriv);
  gamma_deriv.print("Gamma deriv: ");

  // calculate nuclear gradient
  arma::mat grad;
  grad.resize(3, info.atoms.size());
  grad.zeros();
  nuclearGradient(MoleculeAOs, info, grad);
  grad.print("Nuclear gradient: ");


  // calculate X mat
  arma::mat X_mat;
  X_mat.resize(MoleculeAOs.size(), MoleculeAOs.size());
  X_mat.zeros();
  generateXmat(MoleculeAOs, Pa, Pb, X_mat);
  X_mat.print("X matrix: ");

  // geneate y mat 
  arma::mat Y_mat;
  Y_mat.resize(info.num_Atoms, info.num_Atoms);
  Y_mat.zeros();
  generateYmat(MoleculeAOs, info, Pa, Pb, Y_mat);
  Y_mat.print("Y matrix: ");

  // generate elec grad
  arma::mat elec_grad;
  elec_grad.resize(3, info.atoms.size());
  elec_grad.zeros();
  electronGradient(MoleculeAOs, info, S_deriv, gamma_deriv, X_mat, Y_mat, elec_grad);
  elec_grad.print("Electron gradient: ");

  // total gradient
  arma::mat total_grad;
  total_grad.resize(3, info.atoms.size());
  total_grad.zeros();
  total_grad = grad + elec_grad;
  total_grad.print("Gradient: ");

  std::string method_name;
  if (method == "h")
  {
    method_name = "Householder";
  }
  else if (method == "g")
  {
    method_name = "Givens";
  }
  else if (method == "gs")
  {
    method_name = "Gram-Schmidt";
  }
  else if (method == "e")
  {
    method_name = "Armadillo";
  }
  else
  {
    method_name = "Unknown";
  }
  
  cout << "Time for SCF using " << method_name << " method: " << time << " seconds" << endl;

  return EXIT_SUCCESS;
}

