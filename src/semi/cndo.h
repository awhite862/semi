#ifndef CNDO_H
#define CNDO_H

namespace Semi {

class cndo : public semi_method {
private:
    size_t nbasis;
    size_t natoms;
    arma::mat P;
    arma::mat P_prev;
    arma::mat C;
    arma::mat C_prev;
    arma::mat F;

public:
    cndo(/** whatever is necessary **/);

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

#endif // CNDO_H
