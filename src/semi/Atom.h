#ifndef ATOM_H
#define ATOM_H

#include <string>

namespace Semi {
class Atom {
public:
	// xyz coordinates
	double x;
	double y;
	double z;

	//charge or element
	double charge;

	//constructor
	Atom(double xValue, double yValue, double zValue, double chargeValue);

	//get element from charge
	std::string getElement();
};
}//namespace Semi

#endif //ATOM_H
