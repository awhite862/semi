#include "Atom.h"

namespace Semi {

Atom::Atom(double _x, double _y, double _z, double _charge) {
    this->x = _x;
    this->y = _y;
    this->z = _z;
    this->charge = _charge;
}

std::string Atom::getElement() {
    int temp = this->charge;
    switch ((int) temp) {
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

} //namespace Semi
