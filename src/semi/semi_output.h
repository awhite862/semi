#ifndef SEMI_OUTPUT_H
#define SEMI_OUTPUT_H

#include <vector>
#include <armadillo>

namespace Semi {

struct output { 
    double total_energy;
    std::vector<double> orb_energies;
    arma::mat C;
};

} // namespace Semi

#endif // SEMI_OUTPUT_H
