#ifndef BASISSET_H
#define BASISSET_H

#include "GTOFunction.h"
#include "STOFunction.h"
#include "CGTOFunction.h"
#include <vector>

namespace Semi {

template<typename bType>
/** \brief Class that represents a basis set by using an vector of basis functions. **/
class BasisSet {
public:
	std::vector<bType> myBasis; ///!< vector representing basis

public:
	/** \brief Default constructor. **/
	BasisSet() { }

	/** \brief Constructor taking in a vector of basis functions.
	    \param _myBasis vector of basis.
	 **/
	BasisSet(std::vector<bType> _myBasis) {
		myBasis = _myBasis;
	}
};

} // namespace Semi

#endif // BASISSET_H
