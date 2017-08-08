#include "run_cndo.h"
#include <semi/Basis/STOFunction.h>
#include <semi/Integral/Cndo.h>
#include <semi/Huckel/HuckelMethod.h>
#include <cstdlib>
#include <armadillo>
#include <cstdlib>
#include <fstream>
#include "semi/Integral/IntegralEvaluator.h"
#include "semi/Integral/Cndo.h"
#include "semi/Basis/STOFunction.h"
#include "semi/Huckel/HuckelMethod.h"
#include "semi/semi_utils.h"
#include <map>

namespace Semi {

void run_cndo(Molecule &mol, parameters &huckel_params, parameters &cndo_params, output &out) {

   std::vector<CGTOFunction> vbasishuckel;
    for (size_t i = 0; i < mol.myMolecule.size(); i++) {
        std::vector<double> rhuckel(3);
        rhuckel[0] = mol.myMolecule[i].x;
        rhuckel[1] = mol.myMolecule[i].y;
        rhuckel[2] = mol.myMolecule[i].z;
        double charge = mol.myMolecule[i].charge;
        QNumber q1s(1, 0, 0);
        vbasishuckel.push_back(CGTOFunction(q1s, 0, 0, 0, rhuckel, charge));
        if (charge - 0.1 > 2) {
            QNumber q2s(2, 0, 0);
            QNumber q2p(2, 1, 0);
            vbasishuckel.push_back(CGTOFunction(q2s, 0, 0, 0, rhuckel, charge));
            vbasishuckel.push_back(CGTOFunction(q2p, 1, 0, 0, rhuckel, charge));
            vbasishuckel.push_back(CGTOFunction(q2p, 0, 1, 0, rhuckel, charge));
            vbasishuckel.push_back(CGTOFunction(q2p, 0, 0, 1, rhuckel, charge));
        }
    }

    BasisSet<CGTOFunction> bsethuckel(vbasishuckel);
    arma::mat Shuckel;
    calculateOverlapMatrixCGTO(bsethuckel, Shuckel);

    std::vector<STOFunction> vbasis;
    for (size_t i = 0; i < mol.myMolecule.size(); i++) {
        std::vector<double> r(3);
        r[0] = mol.myMolecule[i].x;
        r[1] = mol.myMolecule[i].y;
        r[2] = mol.myMolecule[i].z;
        double charge = mol.myMolecule[i].charge;
        if (charge - 0.1 > 2) {
            QNumber q2s(2, 0, 0);
            QNumber q2px(2, 1, 1);
            QNumber q2py(2, 1, -1);
            QNumber q2pz(2, 1, 0);
            vbasis.push_back(STOFunction(q2s, charge, r[0], r[1], r[2], i));
            vbasis.push_back(STOFunction(q2px, charge, r[0], r[1], r[2], i));
            vbasis.push_back(STOFunction(q2py, charge, r[0], r[1], r[2], i));
            vbasis.push_back(STOFunction(q2pz, charge, r[0], r[1], r[2], i));
        }
        else {
            QNumber q1s(1, 0, 0);
            vbasis.push_back(STOFunction(q1s, charge, r[0], r[1], r[2], i));
        }
    }

    BasisSet<STOFunction> bset(vbasis);
    arma::mat S;
    calculateOverlapMatrix(bset, S);
    double kValue = huckel_params.get_value<double>("k");
    double deltaValue = huckel_params.get_value<double>("delta");
    arma::mat C;
    calculateHuckel(Shuckel, kValue, deltaValue, mol, C);

    C.zeros();
    arma::mat fock;
    SCF(bset, C, S, fock);

    out.C = C;
}

} // namespace Semi


