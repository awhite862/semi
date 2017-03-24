/** \brief Class that represents an atom. **/
#ifndef ATOM_H
#define ATOM_H

#include <string>

namespace Semi {
class Atom {
public:
	/** \brief x value of atom. **/
	double x;

	/** \brief y value of atom. **/
	double y;

	/** \brief z value of atom. **/
	double z;

	/** \brief Charge of an atom. Same as atomic number. **/
	double charge;

	/** \brief Constructor. **/
	Atom(double _x, double _y, double _z, double _charge);

	/** \brief Determines string corresponding to an elements charge. **/
	std::string getElement();
};

} //namespace Semi

#endif //ATOM_H
