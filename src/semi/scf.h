#ifndef SCF_H
#define SCF_H

#include "semi_method.h"
#include <iostream>
#include <sstream>
#include <iomanip>

#include <armadillo>
#include <iostream>
#include <cstdlib>
#include "semi_method.h"
#include "semi/Integral/Cndo.h"
#include "semi/Integral/IntegralEvaluator.h"
#include "semi/Basis/BasisSet.h"
#include "semi/Basis/STOFunction.h"
namespace Semi {

void scf(semi_method &m, int conv, size_t max_iter);

} 

#endif // SCF_H
