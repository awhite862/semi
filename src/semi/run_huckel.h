#ifndef SEMI_RUN_HUCKEL_H
#define SEMI_RUN_HUCKEL_H

#include "Structure/Molecule.h"
#include "parameters.h"
#include "semi_output.h"

namespace Semi {

void run_huckel(Molecule &mol, parameters &huckel_params, output &out);

} // namespace Semi

#endif // SEMI_RUN_HUCKEL_H
