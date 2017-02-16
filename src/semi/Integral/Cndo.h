#ifndef CNDO_H
#define CNDO_H


namespace Semi {
//U_uu
double calculateCoreHamiltonian();

//P_AA
double calculateTotalChargeDensity();

//P_uv
double calculateChargeDensity();

//P_uv
double calculateTotalChargeDensity();

//Gamma_uv
double calculateElectronRepulsion();

//V_AB
double calculateNucleurAttraction();

//B_AB
double calculateBondingParameter();

//S_uv
double calculateOverlapMatrix();

}


#endif //CNDO_H
