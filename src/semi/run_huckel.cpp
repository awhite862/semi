#include "run_huckel.h"
#include <semi/Basis/CGTOFunction.h>
#include <semi/Integral/Cndo.h>
#include <semi/Huckel/HuckelMethod.h>

namespace Semi {

void run_huckel(Molecule &mol, parameters &huckel_params, output &out) {

    std::vector<CGTOFunction> vbasis;
    for (size_t i = 0; i < mol.myMolecule.size(); i++) {
        std::vector<double> r(3);
        r[0] = mol.myMolecule[i].x;
        r[1] = mol.myMolecule[i].y;
        r[2] = mol.myMolecule[i].z;
        double charge = mol.myMolecule[i].charge;
        QNumber q1s(1, 0, 0);
        vbasis.push_back(CGTOFunction(q1s, 0, 0, 0, r, charge));
        if (charge - 0.1 > 2) {
            QNumber q2s(2, 0, 0);
            QNumber q2p(2, 1, 0);
            vbasis.push_back(CGTOFunction(q2s, 0, 0, 0, r, charge));
            vbasis.push_back(CGTOFunction(q2p, 1, 0, 0, r, charge));
            vbasis.push_back(CGTOFunction(q2p, 0, 1, 0, r, charge));
            vbasis.push_back(CGTOFunction(q2p, 0, 0, 1, r, charge));
        }
    }

    BasisSet<CGTOFunction> bset(vbasis);
    arma::mat S;
    calculateOverlapMatrixCGTO(bset, S);
    double kValue = huckel_params.get_value<double>("k");
    double deltaValue = huckel_params.get_value<double>("delta");
    arma::mat C;
    calculateHuckel(S, kValue, deltaValue, mol, C);
    out.C = C;
}

} // namespace Semi
