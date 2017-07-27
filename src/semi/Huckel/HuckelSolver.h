#ifndef HUCKELSOLVER_H
#define HUCKELSOLVER_H
#include "HuckelMethod.h"
#include "semi/Structure/Molecule.h"
#include "semi/semi_utils.h"
#include <armadillo>

namespace Semi {
/** \brief Class that solved Huckel Hamiltonain via diagonalization. **/

/** \brief Solves the Huckel Hamiltonian eigenvalue equation with overlap matrix, HΨ=SΨE. **/
void solveHuckelMatrixWithOverlap(arma::mat Hhuckel, arma::mat Shuckel, arma::mat Smatrix, Molecule m, std::vector<myOrbital> valenceOrbitalData, arma::mat &solMatrix);

/** \brief Solves the Huckel Hamiltonian eigenvalue equation without overlap matrix, HΨ=ΨE. **/
void solveHuckelMatrix(arma::mat Hhuckel, arma::mat Smatrix, Molecule m, std::vector<myOrbital> valenceOrbitalData, arma::mat &solMatrix);



} // namespace Semi

#endif // HuckelSolver.h
