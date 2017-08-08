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
int run_sto_vs_cgto() {
    int i = 0;
    Molecule m;
    std::ifstream fin;
    std::string line;
    std::vector<STOFunction> STOBasis;
    std::vector<CGTOFunction> CGTOBasis;
    fin.open("xyz/o2.txt");
    while (std::getline(fin, line)) {
        std::stringstream linestream(line);
        double elem, x, y, z;
        linestream >> elem >> x >> y >> z;
        std::vector<double> r(3);
        r[0] = x;
        r[1] = y;
        r[2] = z;
        arma::colvec temp;
        temp << x << y << z;
        double charge = elem;
        if (charge - 0.1 > 2) {
            QNumber q2s(2, 0, 0);
            QNumber q2px(2, 1, 1);
            QNumber q2py(2, 1, -1);
            QNumber q2pz(2, 1, 0);
            STOBasis.push_back(STOFunction(q2s, charge, r[0], r[1], r[2], i));
            STOBasis.push_back(STOFunction(q2px, charge, r[0], r[1], r[2], i));
            STOBasis.push_back(STOFunction(q2py, charge, r[0], r[1], r[2], i));
            STOBasis.push_back(STOFunction(q2pz, charge, r[0], r[1], r[2], i));
            CGTOBasis.push_back(CGTOFunction(q2s, 0, 0, 0, temp, elem));
            CGTOBasis.push_back(CGTOFunction(q2px, 1, 0, 0, temp, elem));
            CGTOBasis.push_back(CGTOFunction(q2py, 0, 1, 0, temp, elem));
            CGTOBasis.push_back(CGTOFunction(q2pz, 0, 0, 1, temp, elem));
        }
        else {
            QNumber q1s(1, 0, 0);
            STOBasis.push_back(STOFunction(q1s, charge, r[0], r[1], r[2], i));
            CGTOBasis.push_back(CGTOFunction(q1s, 0, 0, 0, temp, elem));
        }
        i++;
    }
    fin.close();
    BasisSet<STOFunction> STOBasisSet(STOBasis);
    BasisSet<CGTOFunction> CGTOBasisSet(CGTOBasis);
    arma::mat overlapSTO;
    calculateOverlapMatrix(STOBasis, overlapSTO);
    arma::mat overlapCGTO;
    calculateOverlapMatrixCGTO(CGTOBasis, overlapCGTO);
    (round(1000 * overlapSTO) / 1000).print("overlap sto");
    (round(1000 * overlapCGTO) / 1000).print("overap cgto");


    return 0;
}

/** Main method to run tests.
 */
int main() {
    return run_sto_vs_cgto();
}
