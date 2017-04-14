#ifndef STOBASIS_H
#define STOBASIS_H
#include "semi/QNumber.h"

namespace Semi {
/** \brief Class that represents a STO Basis. **/
class STOBasis {
public:
    double zeta; ///!< zeta
    QNumber nlm; ///!< quantum numbers
    double x; ///!< x-position
    double y; ///!< y-position
    double z; ///!< z-position
    double id; ///!< id
    
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
    STOBasis(QNumber _nlm, double _zeta, double _x, double _y, double _z, double id);
};

} //namespace Semi

#endif //STOBASIS_H
