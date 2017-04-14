#ifndef SEMI_UTILS_H
#define SEMI_UTILS_H

#include <string>
#include "semi/Atom.h"

namespace Semi {
/** \brief Determines string corresponding to an elements charge. **/
std::string getElement(double charge);

/** \brief Calculates factorial of n. **/
int factorial(int n);

/** \brief Calculates kronecker delta of i and j. **/
double delta(double i, double j);

} // namespace Semi

#endif // SEMI_UTILS_H
