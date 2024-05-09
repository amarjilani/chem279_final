#pragma once
#ifndef AO_H
#define AO_H
#include <iostream>
#include <armadillo>

struct EmpiricalData {
    double Is_As; // 1/2 * (Is + As)
    double Ip_Ap; // 1/2 * (Ip + Ap)
    int beta; // -β
};
// extern
extern std::map<int, int> valence_map;
extern std::map<int, EmpiricalData> elementData;

double Overlap_onedim(double xa, double xb, double alphaa, double alphab, int la, int lb);

struct MolInfo {
  int num_Atoms;
  int num_elec;
  std::vector<int> atoms;
  int valence_elec;
  int p; 
  int q;

  MolInfo() : num_Atoms(0), num_elec(0), valence_elec(0), p(0), q(0) {}
};

class AO
{
    private:
        arma::vec R0;
        arma::uvec lmn;
        arma::vec alpha;
        arma::vec d_coe;
        arma::vec d_coe_norm;
        int len;
        std::string lable;
        int atomicN;
        int atom_idx;
    public:
        AO(arma::vec &R0_input, arma::vec &alpha_input, arma::vec &d_input, arma::uvec &lmn_input, std::string lable_input, int atomicN, int atom_idx);
        ~AO(){}
        void printinfo();
        arma::uvec get_lmn(){ return lmn;}
        arma::vec get_alpha(){ return alpha;}
        arma::vec get_d_coe(){ return d_coe;}
        arma::vec get_R0(){ return R0;}
        int get_len(){ return len;}
        std::string get_lable(){ return lable;}
        int get_atomicN(){ return atomicN;}
        int get_atom_idx(){ return atom_idx;}
};

MolInfo GenerateAOs(std::vector<AO> &AOs, std::string &fname, arma::mat &H_basis, arma::mat &C_basis, arma::mat &N_basis, arma::mat &O_basis, arma::mat &F_basis);
void PrintAOs(std::vector<AO> &AOs);
double Overlap_onedim(double xa, double xb, double alphaa, double alphab, int la, int lb);
double Overlap_3d(arma::vec &Ra, arma::vec &Rb, double alphaa, double alphab, arma::uvec &lmna, arma::uvec &lmnb);

double Eval_Ov_AOs(AO& sh1, AO& sh2);

void Eval_OV_mat(std::vector<AO> &MoleculeAOs, arma::mat &OV_mat);

#endif // AO_H
