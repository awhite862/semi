#ifndef RUN_SCF_H
#define RUN_SCF_H
#include "parameters.h"
#include "semi_output.h"
#include <semi/Basis/STOFunction.h>
#include <semi/Integral/Cndo.h>
#include <semi/Huckel/HuckelMethod.h>
#include <cstdlib>
#include <armadillo>
#include <cstdlib>
#include <fstream>
#include "semi/Integral/IntegralEvaluator.h"
#include "semi/Integral/Cndo.h"
#include "semi/Basis/STOFunction.h"
#include "semi/Huckel/HuckelMethod.h"
#include "semi/semi_utils.h"
#include <map>


namespace Semi {

void run_scf(Molecule &mol, parameters &huckel_params, parameters &scf_params, parameters &cndo_params, output &out);
} // namespace Semi

#endif // SEMI_RUN_SCF_H
