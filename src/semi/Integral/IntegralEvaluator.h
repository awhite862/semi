#ifndef IntegralEvaluator_H
#define IntegralEvaluator_H

#include "semi/Structure/Molecule.h"
#include "semi/Basis/BasisSet.h"
#include "semi/Basis/CGTOFunction.h"
#include "semi/Basis/GTOFunction.h"
#include "semi/Basis/STOFunction.h"
#include "semi/semi_utils.h"
#include <iostream>
#include <armadillo>

namespace Semi {
/** \brief Class that evaluates overlap integrals. **/

double calculateOverlapSTO(STOFunction a, STOFunction b);

/** \brief calculates overlap integral for STO basis from given parameters. **/
double calculateOverlapSTO(double tau, double rho, double kappa, double rho_alpha, double rho_beta, int *a, int *b);

/** \brief calculates overlap integral for STO basis for different positions and zeta values. **/
double calculateOverlapSTOFull(double tau, double rho, double kappa, double rho_alpha, double rho_beta, int *a, int *b);

/** \brief Calcualtes overlap integral for STO basis given the integrals are centered at same position. **/
double calculateOverlapSTOSamePosition(double tau, double rho, double kappa, double rho_alpha, double rho_beta, int *a, int *b);

/** \brief Calcualtes overlap integral for STO basis given identical zeta values. **/
double calculateOverlapSTOSameZeta(double tau, double rho, double kappa, double rho_alpha, double rho_beta, int *a, int *b);

/** \brief calculates overlap integral for GTO basis from given parameters. **/
double calculateOverlapGTO(GTOFunction a, GTOFunction b);

/** \brief calculates overlap integral for GTO basis from given parameters. **/
double calculateOverlapGTOUnnorm(GTOFunction a, GTOFunction b);

/** \brief calculates overlap integral for CGTO basis from given parameters. **/
double calculateOverlapCGTO(CGTOFunction a, CGTOFunction b);

/** \brief calculates coulomb integral from given parameters. **/
double calculateBasicCoulombIntegral(double zeta_a, double zeta_b, double tau, double rho, double kappa, double rho_alpha, double rho_beta, int *a, int *b);

/** \brief calculates coulomb integral for different positions and zeta values. **/
double calculateBasicCoulombIntegralFull(double zeta_a, double zeta_b, double tau, double rho, double kappa, double rho_alpha, double rho_beta, int *a, int *b);

/** \brief Calcualtes coulomb integral given the integrals are centered at same position. **/
double calculateBasicCoulombIntegralSamePosition(double zeta_a, double zeta_b, double tau, double rho, double kappa, double rho_alpha, double rho_beta, int *a, int *b);

/** \brief Calcualtes coulomb integral given the integrals are centered at same position. **/
double calculateBasicCoulombIntegralSameZeta(double zeta_a, double zeta_b, double tau, double rho, double kappa, double rho_alpha, double rho_beta, int *a, int *b);

/** \brief Calcualtes core valence integral. **/
double calculateCoreValenceInteraction(int *a, int *b);

/** \brief Calcualtes electron repulsion integral. **/
double calculateElectrionRepulsionIntegral(int *a, int *b);

/** \brief Calcualtes basic integral. **/
double calculateBasicIntegral(double zeta_a, double zeta_b, double rho, int *a);

/** \brief Calcualtes rotation matrix to the z axis. **/
arma::mat findRotation(double x1, double y1, double z1, double x2, double y2, double z2);

/** \brief Determines if input for CaluculateOverlap is in wrong order. **/
bool isReversed(int *a, int *b);

} //namespace Semi

#endif //IntegralEvaluator_H
