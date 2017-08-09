#ifndef SCF_H
#define SCF_H

#include "semi_method.h"

namespace Semi {

void scf(semi_method &m, int conv, size_t max_iter);

} 

#endif // SCF_H
