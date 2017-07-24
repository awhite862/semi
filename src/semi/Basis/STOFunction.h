#ifndef STOFUNCTION_H
#define STOFUNCTION_H
#include "semi/Structure/QNumber.h"

namespace Semi {
/** \brief Class that represents a STO Basis. **/
class STOFunction {
public:
    double charge;
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
    STOFunction(QNumber _nlm, double _zeta, double _x, double _y, double _z, double id);
};

} // namespace Semi

#endif // STOFUNCTION_H