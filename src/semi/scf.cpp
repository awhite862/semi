#ifndef SCF_H
#define SCF_H

#include "semi_method.h"

namespace Semi {

void scf(semi_method &m, int conv, size_t max_iter) {

    double thresh = std::pow(10.0,-conv);
    size_t i = 0;
    for (i = 0; i < max_iter; i++) {
        double error = m.get_error();
        double energy = m.get_energy();
        if (std::abs(error) < thresh) break;
        std::cout << energy << "    "  << error << std::endl;
    }
    std::cout << energy << "    "  << error << std::endl;

    if (i == max_iter) {
        std::cout << "Warning: Energy is not converged!!" << std::endl;
    }
    m.final_print();
}

} 

#endif // SCF_H
