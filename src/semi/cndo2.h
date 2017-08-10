#ifndef CNDO2_H
#define CNDO2_H
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

namespace Semi {

class cndo2: public semi_method {
private:
    size_t nbasis;
    size_t natoms;
    size_t nocc;
    arma::mat P;
    arma::mat P_prev;
    arma::mat C;
    arma::mat C_prev;
    arma::mat F;
    arma::mat eigvalues;
    arma::mat S;
    arma::mat occs;
    BasisSet<STOFunction> bset;

public:
    cndo2(size_t nbasis, size_t natoms, arma::mat C, BasisSet<STOFunction> bset);

/** \name Implementation of semi_method **/
///@{
    // return the abs(P - P_prev)
    virtual double get_error();
    // compute and return energy
    virtual double get_energy();
    // Diagonalize F, compute C, P, Build F *** this is better
    // or Build F, diagonalize F, compute C, compute  P
    virtual void take_step();
    // print orbital energies, etc
    virtual void final_print();
///@}
};

} // namespace Semi

#endif // CNDO2_H
