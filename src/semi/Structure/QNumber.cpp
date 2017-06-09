#include "QNumber.h"
namespace Semi {

QNumber::QNumber(int _n, int _l, int _m) {
	n = _n;
	l = _l;
	m = _m;
	std::cout << n << " " <<  l << " " << m << std::endl;
}

} // namespace Semi
