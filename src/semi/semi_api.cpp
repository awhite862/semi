#include "semi_api.h"
#include "run_huckel.h"
#include "run_cndo.h"
#include "run_scf.h"
namespace Semi {

namespace {

void print_start_banner(std::ostream &os) {
    os << " ----------------------------------------------------- " << std::endl
       << "    semi++: Simple semi-empirical quantum chemistry" << std::endl
       << " ----------------------------------------------------- " << std::endl;
    os << "     Authors: Alec White, Yuting Chen" << std::endl << std::endl;
}
void print_end_banner(std::ostream &os) {
    os << " ----------------------------------------------------- " << std::endl
       << " ----------------------------------------------------- " << std::endl;
}

void set_huckel_defaults(parameters &huckel_params) {
    huckel_params.set_default("k", "1.25");
    huckel_params.set_default("delta", "0.5");
}

void set_cndo_defaults(parameters &cndo_params) {
    cndo_params.set_default("guess", "huckel");
    cndo_params.set_default("variant", "2");
}

void set_scf_defaults(parameters &scf_params) {
    scf_params.set_default("maxiter", "100");
    scf_params.set_default("convergence", "5");
    scf_params.set_default("method", "cndo");
}

} //unnamed namespace

void run_semi(input &in, output &out, std::ostream &logger) {

    print_start_banner(logger);

    if (in.ctype == HUCKEL) {
        set_huckel_defaults(in.huckel_params);
        run_huckel(in.mol, in.huckel_params, out);
    }
    if (in.ctype == CNDO) {
        set_huckel_defaults(in.huckel_params);
        set_cndo_defaults(in.cndo_params);
        run_cndo(in.mol, in.huckel_params, in.cndo_params, out);
    }
    if (in.ctype == SCF) {
        set_huckel_defaults(in.huckel_params);
        set_cndo_defaults(in.cndo_params);
        set_scf_defaults(in.scf_params);
        run_scf(in.mol, in.huckel_params, in.cndo_params, in.scf_params, out);
    }

    print_end_banner(logger);
}

} //namespace Semi
