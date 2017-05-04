#ifndef MOLECULE_H
#define MOLECULE_H

#include <vector>
#include "Atom.h"

namespace Semi {
/** \brief Class that represents a molecule by using an vector of atoms. **/
class Molecule {
public:
    std::vector<Atom> myMolecule; ///!< Vector of atoms

public:
    /** \brief Default constructor.
     **/
    Molecule() { }

    /** \brief Constructor. 
        \param _myMolecule Vector of atoms
     **/
    Molecule(std::vector<Atom> _myMolecule);
};

} // namespace Semi

#endif // MOLECULE_H
