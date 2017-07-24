#include "STOFunction.h"

namespace Semi {

STOFunction::STOFunction(QNumber _nlm, double _zeta, double _x, double _y, double _z, double _id) {
	charge = _zeta;
	zeta = _zeta;
	nlm = _nlm;
	x = _x;
	y = _y;
	z = _z;
	id = _id;
	switch ((int)zeta) {
	case 1: zeta = 1.2; break;
	case 2: zeta = 0; break;
	case 3: zeta = 1.3 / 2.0; break;
	case 4: zeta = 1.95 / 2.0; break;
	case 5: zeta = 2.6 / 2.0; break;
	case 6: zeta = 3.25 / 2.0; break;
	case 7: zeta = 3.9 / 2.0; break;
	case 8: zeta = 4.55 / 2.0; break;
	case 9: zeta = 5.2 / 2.0; break;
	default: zeta = 0; break;
	}
}

} // namespace Semi
