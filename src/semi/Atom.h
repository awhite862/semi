#ifndef ATOM_H
#define ATOM_H

#include <string>

namespace Semi {
class Atom {
public:
	// xyz coordinatesdouble _x, double _y, double _z
	double x;
	double y;
	double z;

	//charge or element
	double charge;

	//constructor
	Atom(double _x, double _y, double _z, double _charge);

	//get element from charge
	std::string getElement();
};

} //namespace Semi

#endif //ATOM_H
