#include "BasisSet.h"

namespace Semi {

template <class bType>
BasisSet<bType>::BasisSet(std::vector<bType> _myBasis) {
	myBasis = _myBasis;
}

} // namespace Semi
