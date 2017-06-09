#ifndef GTOFUNCTION_H
#define GTOFUNCTION_H

#include "semi/Structure/QNumber.h"
#include "semi/semi_utils.h"
#include <armadillo>

namespace Semi {
/** \brief Class that represents an GTO Basis Function. **/
class GTOFunction {
public:
    QNumber nlm; ///!< Quantum numbers nlm
    double a; ///!< x power
    double b; ///!< y power
    double c; ///!< z power
    double alpha; ///!< alpha
    arma::colvec r; ///!< col vector representing a vector
    double n; ///!< normalization

public:
    /** \brief Constructor.
        \param _a x power.
        \param _b y power.
        \param _z z power.
        \param _alpha alpha.
        \param _r col vector representing a vector.
     **/
    GTOFunction(QNumber _nlm, double _a, double _b, double _c, double _alpha, arma::colvec _r);
};

} // namespace Semi

#endif // GTOFUNCTION_H
