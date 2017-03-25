#define ARMA_DONT_USE_WRAPPER
#include <iostream>
#include <vector>
#include <armadillo>
#include "semi/Atom.h"
#include "semi/Molecule.h"
#include "semi/Basis.h"
#include "semi/BasisSet.h"
#include "semi/Integral/IntegralEvaluator.h"
#include "semi/Integral/Cndo.h"

using namespace arma;

namespace Semi {
//U_uu
double calculateCoreHamiltonian() {
    return 0;
}

//P_AA
double calculateTotalChargeDensity() {
    return 0;
}

//P_uv
double calculateChargeDensity() {
    return 0;
}

//P_uv
double calculateTotalChargeDensityDiff() {
    return 0;
}

//Gamma_uv
double calculateElectronRepulsion() {
    return 0;
}

//V_AB
double calculateNucleurAttraction() {
    return 0;
}

//B_AB
double calculateBondingParameter() {
    return 0;
}

//S_uv
arma::mat calculateOverlapMatrix(BasisSet a, BasisSet b) {
    mat Smatrix(a.myBasis.size(), b.myBasis.size());
    // for (uint k = 0; k < a.size(); k++) {
    //     for (uint l = 0; l < b.size(); l++) {
    //         r = distance(a[k].x, a[k].y, a[k].z, b[l].x, b[l].y, b[l].z);
    //         zeta_average = 0.5 * (a[k].zeta + b[l].zeta);
    //         tau = (a[k].zeta - b[l].zeta) / (a[k].zeta + b[l].zeta);
    //         rho = 0.5 * (a[k].zeta + b[l].zeta) * r;
    //         kappa = 0.5 * (rho + 1 / rho);
    //         rho_alpha = a[k].zeta * r;
    //         rho_beta = b[l].zeta * r;
    //         int aOrbitalType [3] =  {a[k].n, a[k].l, a[k].m};
    //         int bOrbitalType [3] = {b[l].n, b[l].l, b[l].m};
    //         Smatrix(k, l) = CalculateOverlap(tau, rho, kappa, rho_alpha, rho_beta, aOrbitalType, bOrbitalType);
    //         mat rotationMatrix(3, 3) = findRotation(a[k].x, a[k].y, a[k].z, b[l].x, b[l].y, b[l].z);
    //     }
    // }
    return Smatrix;

}

double distance (double x1, double y1, double z1, double x2, double y2, double z2) {
    return sqrt((pow((x1 - x2), 2) + pow((y1 - y2), 2) + pow((z1 - z2), 2)));
}

} //namespace Semi
