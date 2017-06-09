#ifndef CNDO_H
#define CNDO_H

#include "semi/Integral/IntegralEvaluator.h"

namespace Semi {
//U_uu
double calculateCoreHamiltonian();

//P_AA
arma::mat calculateTotalChargeDensity();

//P_uv
arma::mat calculateChargeDensity();

//Gamma_uv
double calculateElectronRepulsion();

//V_AB
double calculateNucleurAttraction();

//B_AB
double calculateBondingParameter();

//S_uv
arma::mat calculateOverlapMatrix(BasisSet<STOFunction> a);


arma::mat calculateOverlapMatrixGTO(BasisSet<GTOFunction> b);
arma::mat calculateOverlapMatrixCGTO(BasisSet<CGTOFunction> b);

//double distance (double x1, double y1, double z1, double x2, double y2, double z2);

}


#endif //CNDO_H
