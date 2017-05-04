#ifndef SEMI_UTILS_H
#define SEMI_UTILS_H

#include <string>
#include "semi/Atom.h"
#include <cmath>
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

double distance (double x1, double y1, double z1, double x2, double y2, double z2);

} // namespace Semi

#endif // SEMI_UTILS_H
