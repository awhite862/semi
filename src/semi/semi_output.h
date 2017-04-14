#ifndef SEMI_OUTPUT_H
#define SEMI_OUTPUT_H

#include <vector>

namespace Semi {

struct output { 
    double total_energy;
    std::vector<double> orb_energies;
    std::vector<double> C;
};

} // namespace Semi

#endif // SEMI_OUTPUT_H
