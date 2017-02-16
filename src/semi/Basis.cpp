#include "Basis.h"

namespace Semi {
//constructor
Basis::Basis(double zetaValue, int nValue, int lValue, int mValue, double xValue, double yValue, double zValue) {
	this->zeta = zetaValue;
	this->n = nValue;
	this->l = lValue;
	this->m = mValue;
	this->x = xValue;
	this->y = yValue;
	this->z = zValue;
}

//get orbital type
std::string Basis::getOrbitalType() {
	std::stringstream ss;
	ss << this->n;
	std::string s = ss.str();
	switch (this->l) {
	case 0: s += "s";
	case 1: s += "p";
	}
	if (l == 1 && m == 0) {
		s += "z";
	}
	return s;
}
}