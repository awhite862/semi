#define ARMA_DONT_USE_WRAPPER
#include "HuckelMethod.h"
#include "HuckelSolver.h"
#include "HuckelConstructor.h"
#include "HuckelBasis.h"
#include "semi/semi_utils.h"
#include <vector>
namespace Semi {
/*
arma::mat constructCBasis(arma::mat c_full, arma::mat Smatrix) {
    arma::mat y = trans(c_full) * Smatrix * c_full;
    arma::mat c_basis = c_full * (invSqrt(y));
    arma::mat id = trans(c_basis) * Smatrix * c_basis;
    return c_basis;
}*/
void calculateHuckel(arma::mat Smatrix, double kValue, double sigmaValue, Molecule m, arma::mat &solMatrix) {
    arma::mat Hhuckel;
    arma::mat Shuckel;
    std::vector<myOrbital> valenceOrbitalData;

    constructHuckelHamiltonian(Smatrix, kValue, sigmaValue, m, Hhuckel, Shuckel, valenceOrbitalData);
    //solveHuckelMatrixWithOverlap(Hhuckel, Shuckel, Smatrix, m, valenceOrbitalData, solMatrix);
    solveHuckelMatrix(Hhuckel, Smatrix, m, valenceOrbitalData, solMatrix);
}

}