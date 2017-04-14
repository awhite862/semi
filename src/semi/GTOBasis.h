#ifndef GTOBASIS_H
#define GTOBASIS_H

#include "QNumber.h"
#include <armadillo>
#include "semi/semi_utils.h"

namespace Semi {
/** \brief Class that represents an GTO Basis. **/
class GTOBasis {
public:
    double a; ///!< x power
    double b; ///!< y power
    double c; ///!< z power
    double alpha; ///!< alpha
    arma::colvec r; ///!< col vector representing a vector
    int l; ///!< a + b + c
    double n; ///!< normalization

public:
    /** \brief Constructor.
        \param _a x power.
        \param _b y power.
        \param _z z power.
        \param _alpha alpha.
        \param _r col vector representing a vector.
     **/
    GTOBasis(double _a, double _b, double _c, double _alpha, arma::colvec _r);
};

} //namespace Semi

#endif //GTOBASIS_H