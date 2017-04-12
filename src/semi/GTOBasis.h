/** \brief Class that represents a basis. **/
#ifndef GTOBASIS_H
#define GTOBASIS_H

#include "QNumber.h"
#include <armadillo>

namespace Semi {
class GTOBasis {
public:
    QNumber nlm; ///!< zeta
    double N; ///!< n quantum number
    double alpha; ///!< l quantum number
    arma::colvec r; ///!< m quantum number

public:
    /** \brief Constructor.
        \param _zeta zeta.
        \param _n n quantum number.
        \param _l l quantum number.
        \param _m m quantum number.
        \param _x x-position.
        \param _y y-position.
        \param _z z-position.
        \param _id id.
     **/
    GTOBasis(QNumber _nlm, double _N, double _alpha, arma::colvec _r);
};

} //namespace Semi

#endif //GTOBASIS_H
