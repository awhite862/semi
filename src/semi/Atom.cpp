#include "Atom.h"

namespace Semi {

Atom::Atom(double _x, double _y, double _z, double _charge, unsigned _id) {
    x = _x;
    y = _y;
    z = _z;
    charge = _charge;
    id = _id;
}

} //namespace Semi
