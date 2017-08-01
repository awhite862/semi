#ifndef SEMI_RUN_CNDO_H
#define SEMI_RUN_CNDO_H


#include "Structure/Molecule.h"
#include "parameters.h"
#include "semi_output.h"

namespace Semi {

void run_cndo(Molecule &mol, parameters &huckel_params, output &out);

} // namespace Semi

#endif // SEMI_RUN_HUCKEL_H
