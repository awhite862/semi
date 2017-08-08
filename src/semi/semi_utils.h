#ifndef SEMI_UTILS_H
#define SEMI_UTILS_H

#include <string>
#include <cmath>
#include <armadillo>
#include "semi/Basis/BasisSet.h"
#include "semi/Basis/CGTOFunction.h"
#include "semi/Basis/GTOFunction.h"
#include "semi/Basis/STOFunction.h"

namespace Semi {

/** \brief Determines string corresponding to an elements charge. **/
std::string getElement(double charge);

/** \brief Calculates factorial of n. **/
int factorial(int n);

/** \brief Calculates doublefactorial of n. **/
int doubleFactorial(int n);

/** \brief Calculates zeta value for STO. **/
double zetaCalc(double charge);

double delta(double i, double j);

/** \brief Tolerance constant. **/
const double tolerance = 1e-14;

/** \brief Tolerance constant. **/
const double HARTREE_TO_EV = 27.21138602;

/** \brief Calculates distance between points represented by (x1, y1, z1) and (x2, y2, z2). **/
double distance (double x1, double y1, double z1, double x2, double y2, double z2);

/** \brief Calculates the inverse square root of a matrix. **/
void invSqrt(arma::mat A, arma::mat &sol);

/** \brief Calculates number of valence electrons for given atom. **/
int numValence(int num);



/** \brief Fetches Ioniation energy of atom a. **/
double getIE(double a);

/** \brief Fetches electron affinity of atom a. **/
double getEA(double a);

//int numPi(BasisSet<STOFunction> a);

} // namespace Semi

#endif // SEMI_UTILS_H
