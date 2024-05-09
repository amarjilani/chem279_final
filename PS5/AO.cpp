#include "AO.h"
#include <stdexcept>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <sstream>
#include <cassert>
#include <string>
#include <vector>
#include "util.h"
#include <cmath> 
#include <map>

using namespace std;

// valence map
std::map<int, int> valence_map = {{1, 1}, {6, 4}, {7, 5}, {8, 6}, {9, 7}};
std::map<int, EmpiricalData> elementData = {
  {1, {7.176, 0, 9}},
  {6, {14.051, 5.572, 21}},
  {7, {19.316, 7.275, 25}},
  {8, {25.390, 9.111, 31}},
  {9, {32.272, 11.080, 39}}
};


// return number of electrons
MolInfo GenerateAOs(vector<AO> &AOs, string &fname, arma::mat &H_basis, arma::mat &C_basis, arma::mat &N_basis, arma::mat &O_basis, arma::mat &F_basis)
{
  int basislen = H_basis.n_rows;
  assert(C_basis.n_rows == basislen);
  int num_charge;
  int num_Atoms;
  MolInfo molinfo;

  ifstream in(fname, ios::in);
  // cout << fname;

  string line;
  getline(in, line);
  istringstream iss(line);
  if (!(iss >> num_Atoms >> num_charge))
    throw invalid_argument("There is some problem with AO format.");
  int count_atoms = 0;

  while (getline(in, line))
  {
    istringstream iss(line);
    arma::vec R0(3);
    int AtomicN = 0;
    if (!(iss >> AtomicN >> R0[0] >> R0[1] >> R0[2]))
      throw invalid_argument("There is some problem with AO format.");
    
    arma::uvec lmn = {0, 0, 0};
    arma::vec alpha(basislen);
    arma::vec d_coe(basislen);
    // additional basis sets for other atom types 
    if(AtomicN == 1){
      alpha = H_basis.col(0);
      d_coe = H_basis.col(1);
      string lable("H1s");
      AO readedAO(R0, alpha, d_coe, lmn, lable, AtomicN, count_atoms);
      AOs.push_back(readedAO);
    } else if(AtomicN == 6){
        alpha = C_basis.col(0);
        d_coe = C_basis.col(1);
        string lable("C2s");
        AO readedAO(R0, alpha, d_coe, lmn, lable, AtomicN, count_atoms);
        AOs.push_back(readedAO);
        for(size_t j = 0; j < 3; j++){
          d_coe = C_basis.col(2);
          lmn.zeros();
          lmn(j) = 1;
          string lable("C2p");
          AO readedAOp(R0, alpha, d_coe, lmn, lable, AtomicN, count_atoms);
          AOs.push_back(readedAOp);
        }
    } else if (AtomicN == 7){
      alpha = N_basis.col(0);
      d_coe = N_basis.col(1);
      string lable("N2s");
      AO readedAO(R0, alpha, d_coe, lmn, lable, AtomicN, count_atoms);
      AOs.push_back(readedAO);
      for(size_t j = 0; j < 3; j++){
        d_coe = N_basis.col(2);
        lmn.zeros();
        lmn(j) = 1;
        string lable("N2p");
        AO readedAOp(R0, alpha, d_coe, lmn, lable, AtomicN, count_atoms);
        AOs.push_back(readedAOp);
      }
    } else if (AtomicN == 8){
      alpha = O_basis.col(0);
      d_coe = O_basis.col(1);
      string lable("O2s");
      AO readedAO(R0, alpha, d_coe, lmn, lable, AtomicN, count_atoms);
      AOs.push_back(readedAO);
      for(size_t j = 0; j < 3; j++){
        d_coe = O_basis.col(2);
        lmn.zeros();
        lmn(j) = 1;
        string lable("O2p");
        AO readedAOp(R0, alpha, d_coe, lmn, lable, AtomicN, count_atoms);
        AOs.push_back(readedAOp);
      }
    } else if (AtomicN == 9){
      alpha = F_basis.col(0);
      d_coe = F_basis.col(1);
      string lable("F2s");
      AO readedAO(R0, alpha, d_coe, lmn, lable, AtomicN, count_atoms);
      AOs.push_back(readedAO);
      for(size_t j = 0; j < 3; j++){
        d_coe = F_basis.col(2);
        lmn.zeros();
        lmn(j) = 1;
        string lable("F2p");
        AO readedAOp(R0, alpha, d_coe, lmn, lable, AtomicN, count_atoms);
        AOs.push_back(readedAOp);
      }
    } else {
        throw invalid_argument("There are AOs other than HCNOF.");
    }
    // cout << readAO << std::endl;
    count_atoms ++;
    molinfo.atoms.push_back(AtomicN);
    molinfo.valence_elec += valence_map[AtomicN];
  }
  if(count_atoms != num_Atoms){
    throw invalid_argument("Number of AOs are not consistent ");
  }
  molinfo.num_Atoms = num_Atoms;
  molinfo.num_elec = AOs.size() - num_charge;
  if (molinfo.valence_elec % 2 != 0){
    molinfo.p = molinfo.valence_elec / 2 + 1;
    molinfo.q = molinfo.valence_elec / 2;
  } else {
    molinfo.p = molinfo.valence_elec / 2;
    molinfo.q = molinfo.valence_elec / 2;
  }
  in.close();
  return molinfo;
}

// void AO::Reset(double x0_input, double y0_input, double z0_input, double alpha_input, int l_input){
//   R0(0)=x0_input; R0(1)=y0_input; R0(2)=z0_input; alpha=alpha_input; l=l_input;
// }
void PrintAOs(std::vector<AO> &AOs){
  for(auto ao : AOs)
    ao.printinfo();
}


void AO::printinfo(){
  printf("This AO info: %s, R( %1.2f, %1.2f, %1.2f), with angular momentum: %lld %lld %lld\n", lable.c_str(),
    R0(0), R0(1), R0(2), lmn(0), lmn(1), lmn(2));
  d_coe.print("d_coe");
  alpha.print("alpha");
}


double Overlap_onedim(double xa, double xb, double alphaa, double alphab, int la, int lb)
{
  // double x = Combination(3, 1);
  // double y = DoubleFactorial(4);
  double prefactor = exp( -alphaa*alphab*(xa-xb)*(xa-xb)/(alphaa+ alphab))* sqrt(M_PI / (alphaa+ alphab)) ;
  double xP = (alphaa* xa + alphab * xb)/ (alphaa+ alphab);

  double result = 0.0;
  for(int i_index = 0; i_index <= la; i_index++)
    for(int j_index = 0; j_index <= lb; j_index++){
      if((i_index + j_index) % 2 == 1)
        continue;
      double C_part = Combination(la, i_index) * Combination(lb, j_index);
      double DF_part = DoubleFactorial(i_index + j_index - 1);
      double numerator = pow(xP-xa, la -i_index) * pow(xP-xb, lb - j_index);
      // Caution: convert i_index + j_index to float!
      double dominator = pow(2*(alphaa+ alphab), double(i_index + j_index ) / 2.0);
      double temp = C_part * DF_part * numerator / dominator;
      result += temp;
      // printf("%f %f %f %f  %f\n", C_part, DF_part, numerator, dominator, temp);
    }

  result *= prefactor;
  return result;
}

double Overlap_3d(arma::vec &Ra, arma::vec &Rb, double alphaa, double alphab, arma::uvec &lmna, arma::uvec &lmnb){
    double Overlap = Overlap_onedim(Ra(0), Rb(0), alphaa, alphab, lmna(0), lmnb(0)) *
            Overlap_onedim(Ra(1), Rb(1), alphaa, alphab, lmna(1), lmnb(1)) * 
            Overlap_onedim(Ra(2), Rb(2), alphaa, alphab, lmna(2), lmnb(2));
  return Overlap;
}


AO::AO(arma::vec &R0_input, arma::vec &alpha_input,  arma::vec &d_input, arma::uvec &lmn_input, string lable_input, int atomicN, int atom_idx):
R0(R0_input), alpha(alpha_input), d_coe(d_input), lmn(lmn_input), lable(lable_input), atomicN(atomicN), atom_idx(atom_idx){
    assert(R0.n_elem == 3);
    assert(lmn.n_elem == 3);
    len = alpha.n_elem;
    assert(d_coe.n_elem == len);
    for (size_t k = 0; k <len; k++){
      double Overlap_Self = Overlap_3d(R0, R0, alpha(k), alpha(k), lmn, lmn);
      d_coe(k) /= sqrt(Overlap_Self);
    }
}


double Eval_Ov_AOs(AO& sh1, AO& sh2){

  int len = sh1.get_len();
  assert(sh2.get_len() == len);
  arma::vec alphaa = sh1.get_alpha(), alphab = sh2.get_alpha();
  arma::vec Ra = sh1.get_R0(), Rb = sh2.get_R0();
  arma::uvec la = sh1.get_lmn(), lb = sh2.get_lmn();
  arma::vec da = sh1.get_d_coe(), db = sh2.get_d_coe();

  double sum =0.;
  for (size_t k = 0; k <len; k++){
    double alpha_k = alphaa(k);
    for (size_t j = 0; j <len; j++){
      double alpha_j = alphab(j);
      double Overlap = Overlap_3d(Ra, Rb, alpha_k, alpha_j, la, lb);
      // printf("%ld %ld = %1.10f\n", k, j, Overlap);
      sum += da(k) * db(j) * Overlap;
    }
  }
  return sum;
}

void Eval_OV_mat(vector<AO> &MoleculeAOs, arma::mat &OV_mat){
  int dim = MoleculeAOs.size();
  assert(OV_mat.n_rows == dim && OV_mat.n_cols == dim);
  for (size_t k = 0; k <dim; k++){
    for (size_t j = 0; j <= k; j++){
      double OV_1AO = Eval_Ov_AOs(MoleculeAOs[k], MoleculeAOs[j]);
      OV_mat(k,j) = OV_1AO;
      OV_mat(j,k) = OV_1AO;
    }
  }
}

