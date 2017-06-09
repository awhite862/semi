#ifndef HUCKELMETHOD_H
#define HUCKELMETHOD_H
#include "semi/Structure/Molcule.h"
#include <armadillo>
namespace Semi {
/** \brief Class that represents approxiamtes initial guess via Huckel theory. **/
struct voie {
    double atom;
    double oneS;
    double twoS;
    double twoP;
    double threeS;
    double threeP;
};

struct myOrbital {
    Atom atom;
    std::string orbital;
};

/** \brief Huckel approximation of initial guess given overlap matrix, xyz data and parameters.
**/
void calculateHuckel(arma::mat Smatrix, double kValue, double deltaValue, Molecule m, arma::mat &solutionMatrix);

} // namespace Semi

#endif // HUCKELMETHOD_H
