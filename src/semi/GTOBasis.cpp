#include "GTOBasis.h"
#include <math.h>

namespace Semi {
GTOBasis::GTOBasis(QNumber _nlm, double _a, double _b, double _c, double _alpha, arma::colvec _r) {
	nlm = _nlm;
	a = _a;
	b = _b;
	c = _c;
	alpha = _alpha;
	r = _r;
	n = pow(pow(M_PI / (2 * alpha), 3.0 / 2.0) * (doubleFactorial(2.0 * a - 1) * doubleFactorial(2.0 * b - 1) * doubleFactorial(2.0 * c - 1)) / (pow(2, 2 * nlm.l) * pow(alpha, nlm.l)), -1.0 / 2.0);
}

} // namespace Semi
