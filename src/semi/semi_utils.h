#ifndef SEMI_UTILS_H
#define SEMI_UTILS_H

#include <string>
#include "semi/Atom.h"

namespace Semi {
/** \brief Determines string corresponding to an elements charge. **/
std::string getElement(double charge);

/** \brief Calculates factorial of n. **/
int factorial(int n);

/** \brief Calculates doublefactorial of n. **/
int doubleFactorial(int n);

/** \brief Calculates kronecker delta of i and j. **/
double delta(double i, double j);

/** \brief Tolerance constant. **/
const double tolerance = 0.0000000001;

} // namespace Semi

#endif // SEMI_UTILS_H
