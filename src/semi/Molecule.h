#ifndef MOLECULE_H
#define MOLECULE_H

#include <vector>
#include "Atom.h"

namespace Semi {
/** \brief Class that represents a molecule by using an vector of atoms. **/
class Molecule {
public:
    /** \brief vector of basis representing a molecule. **/
    std::vector<Atom> myMolecule;

    /** \brief Constructor. **/
    Molecule(std::vector<Atom> _myMolecule);

};

}//namespace Semi

#endif //MOLECULE_H
