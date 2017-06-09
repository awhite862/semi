#ifndef SEMI_INPUT_H
#define SEMI_INPUT_H

#include "parameters.h"
#include "semi/Structure/Molecule.h"

namespace Semi {

enum calc_type {HUCKEL, CNDO};

struct input {
    Molecule mol; 
    calc_type ctype;
    parameters huckel_params;
    parameters cndo_params;
};

} // namespace Semi

#endif // SEMI_INPUT_H
