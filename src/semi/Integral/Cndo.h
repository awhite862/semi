#ifndef CNDO_H
#define CNDO_H

#include "semi/Integral/IntegralEvaluator.h"

namespace Semi {
//U_uu
double calculateCoreHamiltonian(STOFunction a, STOFunction b);

//P_AA
double calculateTotalChargeDensity(arma::mat density, double id, BasisSet<STOFunction> a) ;

//P_uv
void calculateChargeDensity(arma::mat c_v, arma::mat &density);

//Gamma_uv
double calculateElectronRepulsion(STOFunction a, STOFunction b);

//V_AB
double calculateNucleurAttraction(STOFunction a, STOFunction b);

//B_AB
double calculateBondingParameter(double a, double b) ;

//S_uv
void calculateOverlapMatrix(BasisSet<STOFunction> a, arma::mat &Smatrix);

void SCF(BasisSet<STOFunction> a, arma::mat coefMatrix, arma::mat S, arma::mat &fock);

void calculateOverlapMatrixSTO(BasisSet<STOFunction> b, arma::mat &Smatrix);
void calculateFockMatrix(BasisSet<STOFunction> a, arma::mat coefMatrix, arma::mat S, arma::mat &fock);
void calculateOverlapMatrixGTO(BasisSet<GTOFunction> b, arma::mat &Smatrix);
void calculateOverlapMatrixCGTO(BasisSet<CGTOFunction> b, arma::mat &Smatrix);




//V_AB integral alternative appr
double calculateElectronRepulsionMatagaNishimoto(STOFunction a, STOFunction b);

double calculateElectronRepulsionMatagaNishimotoOffDiag(STOFunction a, STOFunction b) ;
//double distance (double x1, double y1, double z1, double x2, double y2, double z2);

double getIE(double a);
double getEA(double a);
}


#endif //CNDO_H
