#ifndef SEMI_UTILS_H
#define SEMI_UTILS_H

#include <string>
#include <cmath>
#include <armadillo>
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
const double tolerance = 0.0000000001;

/** \brief Calculates distance between points represented by (x1, y1, z1) and (x2, y2, z2). **/
double distance (double x1, double y1, double z1, double x2, double y2, double z2);

/** \brief Calculates the inverse square root of a matrix. **/
void invSqrt(arma::mat A, arma::mat &sol);

} // namespace Semi

#endif // SEMI_UTILS_H
