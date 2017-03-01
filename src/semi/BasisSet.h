#ifndef BASISSET_H
#define BASISSET_H

#include <vector>
#include "Basis.h"

namespace Semi {

class BasisSet {
public:
	std::vector<Basis> myBasis;

	//constructor
	BasisSet(std::vector<Basis> _myBasis);

};

}//namespace Semi
#endif //BASISSET_H
