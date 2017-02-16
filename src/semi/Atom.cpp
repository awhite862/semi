#include "Atom.h"

namespace Semi {
//constructor
Atom::Atom(double xValue, double yValue, double zValue, double chargeValue) {
   this->x = xValue;
   this->y = yValue;
   this->z = zValue;
   this->charge = chargeValue;
}

//get element from charge
std::string Atom::getElement() {
	int temp = this->charge;
	switch((int) temp){
		case 1: return "H";
		case 2: return "He";
		case 3: return "Li";
		case 4: return "Be";
		case 5: return "B";
		case 6: return "C";
		case 7: return "N";
		case 8: return "O";
		case 9: return "F";
		case 10: return "Ne";
	}
	return "null";
}	
}