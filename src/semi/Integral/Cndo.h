#ifndef CNDO_H
#define CNDO_H

#include <armadillo>

namespace Semi {
//U_uu
double calculateCoreHamiltonian();

//P_AA
double calculateTotalChargeDensity();

//P_uv
double calculateChargeDensity();

//P_uv
double calculateTotalChargeDensityDiff();

//Gamma_uv
double calculateElectronRepulsion();

//V_AB
double calculateNucleurAttraction();

//B_AB
double calculateBondingParameter();

//S_uv
arma::mat calculateOverlapMatrix();

double distance (double x1, double y1, double z1, double x2, double y2, double z2);

}


#endif //CNDO_H
