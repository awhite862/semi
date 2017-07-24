/** \brief A test class for semi.*/
#include <armadillo>
#include <cstdlib>
#include <fstream>
#include "semi/Integral/IntegralEvaluator.h"
#include "semi/Integral/Cndo.h"
#include "semi/Basis/CGTOFunction.h"
#include "semi/Huckel/HuckelMethod.h"
#include "semi/semi_utils.h"
#include <map>


using namespace Semi;
/** \brief Test for huckel theory initial guess.*/
int run_huckel_test() {
    arma::mat SMatrix;
    //SMatrix.load("s.txt", arma::raw_ascii);

    int i = 0;
    Molecule m;
    std::ifstream fin;
    std::string line;
    std::vector<CGTOFunction> vbasis;
    fin.open("methane.txt");
    while (std::getline(fin, line)) {
        std::stringstream linestream(line);
        double elem, x, y, z;
        linestream >> elem >> x >> y >> z;
        Atom a(x, y, z, elem, i);
        m.myMolecule.push_back(a);
        double scale = 1.88973;
        std::vector<double> r(3);
        r[0] = x * scale;
        r[1] = y * scale;
        r[2] = z * scale;
        double charge = elem;
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
        i++;
    }
    fin.close();

    //Semi::calculateHuckel(SMatrix, 1, 0.2, Molecule(m.myMolecule), "c_v");

    BasisSet<CGTOFunction> bset(vbasis);
    arma::mat S;
    calculateOverlapMatrixCGTO(bset, S);
    S.print("overlap");
    arma::mat sol;
    Semi::calculateHuckel(S, 1, 0.2, Molecule(m.myMolecule), sol);
    //sol.print();
    return 0;
}

/** Main method to run tests.
 */
int main() {
    return run_huckel_test();
}
