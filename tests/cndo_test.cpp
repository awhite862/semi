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
    std::vector<STOFunction> vbasis;
    fin.open("co.txt");
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

    ///////////////////////////////////////////////////////////////////////////////
    std::vector<CGTOFunction> vbasisprime;
    fin.open("co.txt");
    while (std::getline(fin, line)) {
        std::stringstream linestream(line);
        double elem, x, y, z;
        linestream >> elem >> x >> y >> z;
        Atom a(x, y, z, elem, i);
        double scale = 1;
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
            vbasisprime.push_back(CGTOFunction(q2s, 0, 0, 0, r, charge));
            vbasisprime.push_back(CGTOFunction(q2px, 1, 0, 0, r, charge));
            vbasisprime.push_back(CGTOFunction(q2py, 0, 1, 0, r, charge));
            vbasisprime.push_back(CGTOFunction(q2pz, 0, 0, 1, r, charge));
        }
        else {
            QNumber q1s(1, 0, 0);
            vbasisprime.push_back(CGTOFunction(q1s, 0, 0, 0, r, charge));
        }
        i++;
    }
    fin.close();

    BasisSet<CGTOFunction> bsetprime(vbasisprime);
    arma::mat Sprime;
    calculateOverlapMatrixCGTO(bsetprime, Sprime);
    ///////////////////////////////////////////////////////////////////////
    BasisSet<STOFunction> bset(vbasis);
    arma::mat S;
    calculateOverlapMatrixSTO(bset, S);
    //arma::mat Srot = calculateOverlapMatrix(bset);


    arma::mat sol;
    Semi::calculateHuckel(S, 1, 0.2, Molecule(m.myMolecule), sol);
    sol.print("sol");

    arma::mat fock;

    // Sprime(0, 4) = 0.4873;
    // Sprime(4, 0) = 0.4873;
    // Sprime(3, 4) = 0.5172;
    // Sprime(4, 3) = 0.5172;
    //Sprime.eye(4,4);


    //S.load("bh2.txt", arma::raw_ascii);
    S.print("Overlap matrix");
    sol.zeros();
    sol.print("coeffs");
    SCF(bset, sol, S, fock);





    arma::mat colbyeigvec;
    colbyeigvec.load("colby.txt", arma::raw_ascii);
    //colbyeigvec.print("colby");

    arma::mat colbyeigval;
    colbyeigval.load("colby.txt", arma::raw_ascii);
    //colbyeigval.print("colby");

    arma::mat temp;
    temp = colbyeigvec * diagmat(colbyeigval) * inv(colbyeigvec);
    temp.print("fock");
    //Sprime.print("Overlap matrix cgto");

    S.print("Overlap matrix sto");

    return 0;
}

/** Main method to run tests.
 */
int main() {
    return run_huckel_test();
}
