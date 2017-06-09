#include "Molecule.h"
#include "semi/semi_utils.h"
namespace Semi {

Molecule::Molecule(std::vector<Atom> _myMolecule) {
	myMolecule = _myMolecule;
}

double Molecule::getNuclearEnergy() {
	double sum = 0;
	for (int k = 0; k < myMolecule.size(); k++) {
		for (int j = 0; j < myMolecule.size(); j++) {
			double dist = distance(myMolecule[k].x, myMolecule[k].y, myMolecule[k].z, myMolecule[j].x, myMolecule[j].y, myMolecule[j].z);
			if (dist > tolerance) {
				sum += myMolecule[k].charge * myMolecule[j].charge / dist;
			}
		}
	}
}

} // namespace Semi
