#ifndef BASISSET_H
#define BASISSET_H

#include <vector>
#include "GTOBasis.h"
#include "STOBasis.h"
#include "CGTOBasis.h"

namespace Semi {
template<typename bType>
/** \brief Class that represents a basis set by using an vector of basis. **/
class BasisSet {
public:
	std::vector<bType> myBasis; ///!< vector representing basis

public:
	/** \brief Constructor.
	    \param _myBasis vector of basis.
	 **/
	BasisSet(std::vector<bType> _myBasis) {
            myBasis = _myBasis;
        }
};

} //namespace Semi

#endif //BASISSET_H
