#ifndef SEMI_INPUT_H
#define SEMI_INPUT_H

#include "parameters.h"
#include "Molecule.h"

namespace Semi {

enum calc_type {HUCKEL, CNDO};

struct input {
    Molecule mol; 
    parameters huckel_params;
    parameters cndo_params;
};

} // namespace Semi

#endif // SEMI_INPUT_H
