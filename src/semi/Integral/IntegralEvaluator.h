/** \brief Class that evaluates overlap integrals. **/
#ifndef IntegralEvaluator_H
#define IntegralEvaluator_H

#include <iostream>
#include <vector>
#include <armadillo>
#include <semi/Atom.h>
#include <semi/Molecule.h>
#include <semi/Basis.h>
#include <semi/BasisSet.h>
#include "IntegralEvaluator.h"

namespace Semi {
/** \brief Calculates overlap integral from given parameters.
 *  Higher level method that determines which formula to use to calculate overlap integral.
 **/
double CalculateOverlap(double tau, double rho, double kappa, double rho_alpha, double rho_beta, int *a, int *b);

/** \brief Calculates overlap integral for different positions and zeta values. **/
double CalculateOverlapFull(double tau, double rho, double kappa, double rho_alpha, double rho_beta, int *a, int *b);

/** \brief Calcualtes overlap integral given the integrals are centered at same position. **/
double CalculateOverlapSamePosition(double tau, double rho, double kappa, double rho_alpha, double rho_beta, int *a, int *b);

/** \brief Calcualtes overlap integral given identical zeta values. **/
double CalculateOverlapSameZeta(double tau, double rho, double kappa, double rho_alpha, double rho_beta, int *a, int *b);

/** \brief Calcualtes rotation matrix to the z axis. **/
arma::mat findRotation(double x1, double y1, double z1, double x2, double y2, double z2);

/** \brief Determines if input for CaluculateOverlap is in wrong order. **/
bool isReversed(int *a, int *b);

/** \brief Tolerance constant. **/
const double tolerance = 0.1;

/**	\brief Calculates factorial given an integer recursively. **/
int factorial(int n);

} //namespace Semi

#endif //IntegralEvaluator_H
