#ifndef BASISSET_H
#define BASISSET_H

#include <vector>
#include "Basis.h"

namespace Semi {
/** \brief Class that represents a basis set by using an vector of basis. **/
class BasisSet {
public:
    /** \brief vector of basis representing a basis set. **/
    std::vector<Basis> myBasis;

    /** \brief Constructor. **/
    BasisSet(std::vector<Basis> _myBasis);

};

}//namespace Semi

#endif //BASISSET_H
