#include "semi_api.h"

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
} // unnamed namespace

void run_semi(input in, output out, std::ostream &logger) {

    print_start_banner(logger);

    print_end_banner(logger);
}

} // namespace Semi
