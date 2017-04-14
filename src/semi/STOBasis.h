#ifndef STOBASIS_H
#define STOBASIS_H

namespace Semi {
/** \brief Class that represents a STO Basis. **/
class STOBasis {
public:
    double zeta; ///!< zeta
    int n; ///!< n quantum number
    int l; ///!< l quantum number
    int m; ///!< m quantum number
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
    STOBasis(double _zeta, int _n, int _l, int _m, double _x, double _y, double _z, double id);
};

} //namespace Semi

#endif //STOBASIS_H
