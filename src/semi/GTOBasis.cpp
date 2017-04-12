#include "GTOBasis.h"

namespace Semi {

GTOBasis::GTOBasis(QNumber _nlm, double _N, double _alpha, arma::colvec _r) {
	nlm = _nlm;
	N = _N;
	alpha = _alpha;
	r = _r;
}

} //namespace Semi
