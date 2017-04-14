#include "QNumber.h"

namespace Semi {
QNumber::QNumber(int _n, int _l, int _m) {
	n = _n;
	l = _l;
	m = _m;
}

QNumber::QNumber() {
	n = 0;
	l = 0;
	m = 0;
}

} //namespace Semi
