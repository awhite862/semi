#include <fstream>
#include <io/io_input_section.h>
#include <cstdlib> // EXIT_FAILURE failed compilation
#include "../src/semi/Atom.h"
#include <armadillo>

namespace Semi {
/** \test basic compilation. **/
int run_compilation_test() {
    Atom* a = new Atom(1,1,1,1);
    return 0;
}


int main() {
    return run_compilation_test() | 0;
}

}