#include "cndo2.h"
#include <armadillo>
#include <iostream>
#include <cstdlib>
#include "semi_method.h"
#include "semi/Integral/Cndo.h"
#include "semi/Integral/IntegralEvaluator.h"
#include "semi/Basis/BasisSet.h"
#include "semi/Basis/STOFunction.h"

namespace Semi {


cndo2::cndo2(size_t _nbasis, size_t _natoms, arma::mat _C, BasisSet<STOFunction>  _bset) {
    nbasis = _nbasis;
    natoms = _natoms;
    C = _C;
    bset = _bset;

    //previous to zeros
    P_prev.set_size(C.n_rows, C.n_rows);
    P_prev.zeros();
    C_prev.set_size(C.n_rows, C.n_rows);
    C_prev.zeros();

    //S
    calculateOverlapMatrix(bset, S);

    unsigned int numpi = 0;
    std::vector<int> ids;
    for (int k = 0; k < bset.myBasis.size(); k++) {
        if (std::find(ids.begin(), ids.end(), bset.myBasis[k].id) == ids.end()) {
            ids.push_back(bset.myBasis[k].id);
            numpi += numValence(bset.myBasis[k].charge);
        }
    }
    nocc = numpi / 2;
    occs = C.cols(0, nocc - 1);
    //P
    calculateChargeDensity(occs, P);
    //Initial guess
    calculateFockMatrix3(bset, occs, S, F);
    eigvalues = C;
}

// return the abs(P - P_prev)
double cndo2::get_error() {
    return norm(P - P_prev);
}
// compute and return energy
double cndo2::get_energy() {
    double energy = 0;
    for (int i = 0; i < nocc; i++) {
        energy += 2.0 * eigvalues(i);
    }
    return energy;
}
// Diagonalize F, compute C, P, Build F *** this is better
// or Build F, diagonalize F, compute C, compute  P
void cndo2::take_step() {
    //update
    P_prev = P;
    C_prev = C;

    //diagonalize
    arma::mat eigvec;
    arma::vec eigval;
    eig_sym(eigval, eigvec, F);
    eigvalues = eigval;
    //compute C
    C = eigvec;
    occs = C.cols(0, nocc - 1);

    //compute P
    calculateChargeDensity(occs, P);

    //Build F
    calculateFockMatrix3(bset, occs, S, F);

}
// print orbital energies, etc
void cndo2::final_print() {

    double energy = 0;
    std::cout << "Orbital Energies" << std::endl;
    for (int i = 0; i < eigvalues.size(); i++) {
        if (eigvalues(i) > 0)
            std::cout <<  std::fixed << "Energy = " << "+" << eigvalues(i) << std::endl;
        else
            std::cout <<  std::fixed << "Energy = " << "-" << (-1.00 * eigvalues(i)) << std::endl;
    }
    for (int i = 0; i < nocc; i++) {
        energy += 2.0 * eigvalues(i);
    }
    std::cout << "Total Energy: " << energy << std::endl;
}
///@}

} //namespace Semi
