#ifndef IntegralEvaluator_H
#define IntegralEvaluator_H

#include "semi/Structure/Molecule.h"
#include "semi/Basis/BasisSet.h"
#include "semi/Basis/CGTOFunction.h"
#include "semi/Basis/GTOFunction.h"
#include "semi/Basis/STOFunction.h"
#include "semi/semi_utils.h"
#include "semi/Structure/QNumber.h"
#include <iostream>
#include <armadillo>

namespace Semi {
/** \brief Class that evaluates overlap integrals. **/

/** \brief Calculates overlap matrix with rotation. **/
void calculateOverlapMatrix(BasisSet<STOFunction> a, arma::mat &Smatrix);

/** \brief Calculates overlap integral for between two CGTO basis functions. **/
double calculateOverlapCGTO(CGTOFunction a, CGTOFunction b);

/** \brief Calculates overlap integral for between two GTO basis functions. **/
double calculateOverlapGTO(GTOFunction a, GTOFunction b);

/** \brief Calculates the unnormalized overlap integral for between two BGO basis functions. **/
double calculateOverlapGTOUnnorm(GTOFunction a, GTOFunction b);

/** \brief Calculates overlap integral for between two STO basis functions. **/
double calculateOverlapSTO(STOFunction a, STOFunction b);

/** \brief Calculates overlap integral for STO basis from given parameters. **/
double calculateOverlapSTO(double tau, double rho, double kappa, double rho_alpha, double rho_beta, QNumber aQNumber, QNumber bQNumber);

/** \brief Calculates overlap integral for STO basis for different positions and zeta values. **/
double calculateOverlapSTOFull(double tau, double rho, double kappa, double rho_alpha, double rho_beta, QNumber aQNumber, QNumber bQNumber);

/** \brief Calcualtes overlap integral for STO basis given identical zeta values. **/
double calculateOverlapSTOSameZeta(double tau, double rho, double kappa, double rho_alpha, double rho_beta, QNumber aQNumber, QNumber bQNumber);

/** \brief Calcualtes overlap integral for STO basis given the integrals are centered at same position. **/
double calculateOverlapSTOSamePosition(double tau, double rho, double kappa, double rho_alpha, double rho_beta, QNumber aQNumber, QNumber bQNumber);

/** \brief Calculates coulomb integral from given parameters. **/
double calculateBasicCoulombIntegral(double zeta, double tau, double rho, double kappa, double rho_alpha, double rho_beta, QNumber aQNumber, QNumber bQNumber);

/** \brief Calculates coulomb integral for different positions and zeta values. **/
double calculateBasicCoulombIntegralFull(double zeta, double tau, double rho, double kappa, double rho_alpha, double rho_beta, QNumber aQNumber, QNumber bQNumber);

/** \brief Calcualtes coulomb integral given the integrals are centered at same position. **/
double calculateBasicCoulombIntegralSameZeta(double zeta, double tau, double rho, double kappa, double rho_alpha, double rho_beta, QNumber aQNumber, QNumber bQNumber);

/** \brief Calcualtes coulomb integral given the integrals are centered at same position. **/
double calculateBasicCoulombIntegralSamePosition(double zeta, double tau, double rho, double kappa, double rho_alpha, double rho_beta, QNumber aQNumber, QNumber bQNumber);

/** \brief Calcualtes rotation matrix to the z axis	to correct overlap integrals. **/
void findRotation(double x1, double y1, double z1, double x2, double y2, double z2, arma::mat &rotation);

/** \brief Determines if input for CaluculateOverlap is in wrong order. **/
bool isReversed(QNumber aQNumber, QNumber bQNumber);

} //namespace Semi

#endif //IntegralEvaluator_H
