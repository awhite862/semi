#include "semi_api.h"
#include "run_huckel.h"
#include "run_cndo.h"
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
} //unnamed namespace

void run_semi(input &in, output &out, std::ostream &logger) {

    print_start_banner(logger);
    
    if (in.ctype == HUCKEL) {
        run_huckel(in.mol, in.huckel_params, out);
    }
    if (in.ctype == CNDO) {
        run_cndo(in.mol, in.huckel_params, out);
    }

    print_end_banner(logger);
}

} //namespace Semi
