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

    std::vector<STOFunction> valenceBasisSet;
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
            valenceBasisSet.push_back(STOFunction(q2s, charge, r[0], r[1], r[2], i));
            valenceBasisSet.push_back(STOFunction(q2px, charge, r[0], r[1], r[2], i));
            valenceBasisSet.push_back(STOFunction(q2py, charge, r[0], r[1], r[2], i));
            valenceBasisSet.push_back(STOFunction(q2pz, charge, r[0], r[1], r[2], i));
        }
        else {
            QNumber q1s(1, 0, 0);
            valenceBasisSet.push_back(STOFunction(q1s, charge, r[0], r[1], r[2], i));
        }
    }
    BasisSet<STOFunction> bset(valenceBasisSet);

    arma::mat overlapMatrix;
    calculateOverlapMatrix(bset, overlapMatrix);

    double kValue = huckel_params.get_value<double>("k");
    double deltaValue = huckel_params.get_value<double>("delta");

    cndo_params.set_default("guess","huckel");
    cndo_params.set_default("max_iter","100");
    cndo_params.set_default("convergence","5");
    cndo_params.set_default("variant","2");

    double guessType = cndo_params.get_value<double>("guess");
    double maxIterations = cndo_params.get_value<double>("max_iter");
    double convergence = cndo_params.get_value<double>("convergence");
    double cndoVariant = cndo_params.get_value<double>("variant");

    arma::mat C;
    //calculateHuckel(Shuckel, kValue, deltaValue, mol, C);
    C.zeros();

    arma::mat fock;
    SCF(bset, C, overlapMatrix, maxIterations, convergence, fock);

    out.C = C;
}

} // namespace Semi


