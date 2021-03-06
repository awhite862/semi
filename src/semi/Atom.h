#ifndef ATOM_H
#define ATOM_H

#include <string>

namespace Semi {
/** \brief Class that represents an atom. **/
class Atom {
public:
    double x; ///!< x-position
    double y; ///!< y-position
    double z; ///!< z-position
    double charge; ///!< charge of the atom

public:
    /** \brief Default constructor
     **/
    Atom() : x(0), y(0), z(0), charge(0.0) { }

    /** \brief Constructor. 
        \param _x x-position.
        \param _y x-position.
        \param _z x-position.
        \param _charge charge.
     **/
    Atom(double _x, double _y, double _z, double _charge);

    /** \brief Determines string corresponding to an elements charge. **/
    std::string getElement();
};

} //namespace Semi

#endif //ATOM_H
