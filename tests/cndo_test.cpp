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

    int i = 0;
    Molecule m;
    std::ifstream fin;
    std::string line;
    std::vector<STOFunction> vbasis;
    fin.open("xyz/o2.txt");
    while (std::getline(fin, line)) {
        std::stringstream linestream(line);
        double elem, x, y, z;
        linestream >> elem >> x >> y >> z;
        Atom a(x, y, z, elem, i);
        m.myMolecule.push_back(a);
        double scale = 1.;
        std::vector<double> r(3);
        r[0] = x * scale;
        r[1] = y * scale;
        r[2] = z * scale;
        double charge = elem;
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
        i++;
    }
    fin.close();
    BasisSet<STOFunction> bset(vbasis);
    arma::mat S;
    calculateOverlapMatrix(bset, S);
    arma::mat sol;
    Semi::calculateHuckel(S, 1, 0.2, Molecule(m.myMolecule), sol);
    arma::mat fock;
    SCF(bset, sol, S, 50, 9, fock);

    std::cout << "------------------------------------------------------------------------------------" << std::endl;
    calculateOverlapMatrix(bset, S);
    (round(1000 * S) / 1000).print("Overlap matrix sto");
    arma::mat eigvec;
    arma::vec eigval;
    eig_sym(eigval, eigvec, S);
    eigval.print("eigval");

    return 0;
}

/** Main method to run tests.
 */
int main() {
    return run_huckel_test();
}
