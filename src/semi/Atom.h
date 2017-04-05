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
    double id; ///!< id for the atom

public:
    /** \brief Constructor. 
        \param _x x-position.
        \param _y y-position.
        \param _z z-position.
        \param _charge charge.
     **/
    Atom(double _x, double _y, double _z, double _charge, double _id);

    /** \brief Determines string corresponding to an elements charge. **/
    std::string getElement();
};

} //namespace Semi

#endif //ATOM_H
