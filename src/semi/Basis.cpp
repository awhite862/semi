#include "Basis.h"

namespace Semi {

//constructor
Basis::Basis(double _zeta, int _n, int _l, int _m, double _x, double _y, double _z) {
	this->zeta = _zeta;
	this->n = _n;
	this->l = _l;
	this->m = _m;
	this->x = _x;
	this->y = _y;
	this->z = _z;
}

} //namespace Semi