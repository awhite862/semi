#ifndef SCF_H
#define SCF_H

#include "semi_method.h"
#include <iostream>
#include <sstream>
#include <iomanip>

#include <armadillo>
#include <iostream>
#include <cstdlib>
#include "semi_method.h"
#include "semi/Integral/Cndo.h"
#include "semi/Integral/IntegralEvaluator.h"
#include "semi/Basis/BasisSet.h"
#include "semi/Basis/STOFunction.h"
#include "semi_method.h"

namespace Semi {

void scf(semi_method &m, int conv, size_t max_iter) {

    std::cout << " ----------------------------------------------------- " << std::endl;
    std::cout << std::setw(10) << "Iteration" << std::setw(20) << "Total Energy" << std::setw(15) << "Error" << std::setw(15) << std::endl;
    std::cout << " ----------------------------------------------------- " << std::endl;
    double thresh = std::pow(10.0, -conv);
    size_t i = 0;
    double energy = 0, error = 0 ;
    for (i = 0; i < max_iter; i++) {
        error = m.get_error();
        energy = m.get_energy();
        if (std::abs(error) < thresh) {
            std::cout << std::setw(10) << std::fixed << i << std::setw(20) << std::setprecision(5) << energy << std::setw(20) << std::scientific << std::setprecision(5) << error << std::setw(10) << "     Convergence criteria met" << std::endl;
            break;
        }
        std::cout << std::setw(10) << std::fixed << i << std::setw(20) << std::setprecision(5) << energy << std::setw(20) << std::scientific << std::setprecision(5) << error << std::setw(10) << std::endl;
        m.take_step();
    }

    if (i == max_iter) {
        std::cout << "Warning: Energy is not converged!!" << std::endl;
    }
    m.final_print();
}

}

#endif // SCF_H
