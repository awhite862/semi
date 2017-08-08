#ifndef CNDO_H
#define CNDO_H

#include "semi/Integral/IntegralEvaluator.h"

namespace Semi {

void printEnergy(arma::mat eigval, double num_orbitals, double variant);
/** \brief Iteratively . **/
void SCF(BasisSet<STOFunction> a, arma::mat coefMatrix, arma::mat S, double maxIterations, double convergence, arma::mat &fock);


///////////////////////////////////////////////////////////////////////////////////////////////////

/** \brief Constructs fock matrix for a basisset a given coefficients and overlap matricies for CNDO(1). F_uv. **/
void calculateFockMatrix(BasisSet<STOFunction> a, arma::mat coefMatrix, arma::mat S, arma::mat &fock);

/** \brief Calculates core Hamiltonian for atom a, for CNDO(I). U_uu **/
double calculateCoreHamiltonian(STOFunction a, STOFunction b);

std::vector<int> calculateNumberValenceElectrons(double charge);
double calculateIonizationPotential(double charge, double l);

/** \brief Calculates charge desnity of atom represented by id in density matrix density. P_AA. **/
double calculateTotalChargeDensity(arma::mat density, double id, BasisSet<STOFunction> a) ;

/** \brief Calculates charge density from coefficient matrix c_v. P_uv. **/
void calculateChargeDensity(arma::mat c_v, arma::mat &density);

/** \brief Calculates electron repulsion between atoms a and b. Gamma_uv. **/
double calculateElectronRepulsion(STOFunction a, STOFunction b);

/** \brief Calculates nuclear attraction between atoms a and b for CNDO(I). V_AB. **/
double calculateNucleurAttraction(STOFunction a, STOFunction b);

/** \brief Calculates bonding parameter between atoms a and b. B_ab. **/
double calculateBondingParameter(double a, double b);

//helper method for B_AB
double getBondingParameter(double a);

/** \brief Constructs fock matrix for a basis set a given coefficients and overlap matricies for CNDO(2). F_uv. **/
void calculateFockMatrix3(BasisSet<STOFunction> a, arma::mat coefMatrix, arma::mat S, arma::mat & fock);

/** \brief Calculates atomic data term for CNDO(2). **/
double calculateAtomicData(double charge, double l);

/** \brief Calculates electron repulsion between atoms a and b using Mataga Nishimoto Approximation. Gamma_uv. **/
double calculateElectronRepulsionMatagaNishimoto(STOFunction a, STOFunction b);

/** \brief Calculates electron repulsion for off diagonal elements between atoms a and b using Mataga Nishimoto Approximation. Gamma_uv. **/
double calculateElectronRepulsionMatagaNishimotoOffDiag(STOFunction a, STOFunction b) ;







/** \brief Calculates overlap matrix for a STO basis set. **/
//void calculateOverlapMatrixSTO(BasisSet<STOFunction> b, arma::mat &Smatrix);


/** \brief Calculates overlap matrix for a CGTO basis set. **/
void calculateOverlapMatrixCGTO(BasisSet<CGTOFunction> b, arma::mat &Smatrix);

}


#endif //CNDO_H
