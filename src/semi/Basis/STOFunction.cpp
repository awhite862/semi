#include "STOFunction.h"

namespace Semi {
	
STOFunction::STOFunction(QNumber _nlm, double _zeta, double _x, double _y, double _z, double _id) {
	zeta = _zeta;
	nlm = _nlm;
	x = _x;
	y = _y;
	z = _z;
	id = _id;
}

} // namespace Semi
