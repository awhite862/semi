#include "run_huckel.h"
#include "CGTOBasis.h"

namespace Semi {

void run_huckel(Molecule &mol, parameters &huckel_params, output &out) {
    
    std::vector<CGTOBasis> vbasis;
    for (size_t i = 0; i < mol.myMolecule.size(); i++) {
        std::vector<double> r(3);
        r[0] = mol.myMolecule[i].x;
        r[1] = mol.myMolecule[i].y;
        r[2] = mol.myMolecule[i].z;
        double charge = mol.myMolecule[i].charge;
        if (charge + 0.1 < 3) { // add 1S
            QNumber q1s(1,0,0);
            vbasis.push_back(CGTOBasis(q1s, 0, 0, 0, r, charge));
        }
        else { // add 2s 2p
            QNumber q2s(2,0,0);
            QNumber q2p(2,1,0);
            vbasis.push_back(CGTOBasis(q2s, 0, 0, 0, r, charge));
            vbasis.push_back(CGTOBasis(q2p, 1, 0, 0, r, charge));
            vbasis.push_back(CGTOBasis(q2p, 0, 1, 0, r, charge));
            vbasis.push_back(CGTOBasis(q2p, 0, 0, 1, r, charge));
        }
    }

}

} // namespace Semi
