#include "GTOBasis.h"
#include <math.h>

namespace Semi {
GTOBasis::GTOBasis(double _a, double _b, double _c, double _alpha, arma::colvec _r) {
	a = _a;
	b = _b;
	c = _c;
	alpha = _alpha;
	r = _r;
	l = a + b + c;
	n = M_PI / (2 * alpha) * pow((factorial(factorial(2.0 * a - 1)) * factorial(factorial(2.0 * b - 1)) * factorial(factorial(2.0 * c - 1))) / (pow(2, 2 * l) * pow(alpha, l)), -1.0 / 2.0);
}

} //namespace Semi
