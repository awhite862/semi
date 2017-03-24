/** \brief A test class for semi.*/
#include <fstream>
#include <io/io_input_section.h>
#include <cstdlib>
#include <semi/Atom.h>
#include <semi/Basis.h>
#include <semi/Molecule.h>
#include <semi/BasisSet.h>
#include <semi/Integral/IntegralEvaluator.h>
#include <semi/Integral/Cndo.h>
#include <semi/Huckel/huckel.h>
//#include <armadillo>
using namespace Semi;
/** \brief Basic compilation test for semi classes.*/
int run_compilation_test() {
	/*Semi::Atom * a = new Semi::Atom(1, 1, 1, 1);
	Semi::Basis * b = new Semi::Basis(1, 1, 1, 1, 1, 1, 1);
	Semi::Molecule * Molecule = new Semi::Molecule(std::vector<Semi::Atom> (4, *a));
	Semi::BasisSet * BasisSet = new Semi::BasisSet(std::vector<Semi::Basis> (4, *b));*/

	return 0;
}

/** \brief Basic Overlap Integral test.
 *
 *  Overlap of 2 1s orbitals, centered at (0, 0, 0), (1, 0, 0)
 *	1st test zeta values 1, 1, evaluates to 0.214596340683341354 un-normalized
 *	2nd test zeta values 0.2, 1, evaluates to 0.121087705563123644 un-normalized
 */
int run_sto_test() {
	///Unnormalized values
	double actualFirst = 0.214596340683341354;
	double actualSecond = 1.121087705563123644;
	double tol = 0.0000000001;

	///paramaters for the first overlap integral
	double r = 1;
	double z1 = 1;
	double z2 = 1;
	double tau = (z1 - z2) / (z1 + z2);
	double rho = 0.5 * (z1 + z2) * r;
	double kappa = 0.5 * (tau + 1 / tau);
	double rho_alpha = z1 * r;
	double rho_beta = z2 * r;
	int a[3];
	a[0] = 1;
	a[1] = 0;
	a[2] = 0;
	int b[3];
	b[0] = 1;
	b[1] = 0;
	b[2] = 0;
	double first = CalculateOverlap(tau, rho, kappa, rho_alpha, rho_beta, a, b);
	double normFirst = pow(2 * z1, 1 + 0.5) * pow(Semi::factorial(2 * 1), -0.5) *  pow(2 * z2, 1 + 0.5) * pow(Semi::factorial(2 * 1), -0.5) * actualFirst;

	///parameters for second overlap integral
	z1 = 0.2;
	z2 = 1;
	tau = (z1 - z2) / (z1 + z2);
	rho = 0.5 * (z1 + z2) * r;
	kappa = 0.5 * (tau + 1.0 / tau);
	rho_alpha = z1 * r;
	rho_beta = z2 * r;
	double second = Semi::CalculateOverlap(tau, rho, kappa, rho_alpha, rho_beta, a, b);
	double normSecond = pow(2 * z1, 1 + 0.5) * pow(Semi::factorial(2 * 1), -0.5) *  pow(2 * z2, 1 + 0.5) * pow(Semi::factorial(2 * 1), -0.5) * actualSecond;

	normFirst -= first;
	normSecond -= second;
	return !(abs(normFirst) < tol || abs(normSecond) < tol);
}

/** \brief Test for rotation matrix in overlap integrals.*/
int run_rotation_test() {
	double tol = 0.0000000001;

	arma::vec vec1(3);
	vec1(0) = 5;
	vec1(1) = 5;
	vec1(2) = 5;
	arma::mat rotated1 = Semi::findRotation(5, 5, 5, 7, 7, 7);

	arma::vec vec2(3);
	vec2(0) = 5;
	vec2(1) = 5;
	vec2(2) = 5;
	arma::mat rotated2 = Semi::findRotation(5, 5, 5, 3, 3, 3);

	arma::vec vec3(3);
	vec3(0) = 1;
	vec3(1) = 1;
	vec3(2) = 1;
	arma::mat rotated3 = Semi::findRotation(1, 1, 1, 3, 4, 5);

	arma::vec vec4(3);
	vec4(0) = 1;
	vec4(1) = -2;
	vec4(2) = -3;
	arma::mat rotated4 = Semi::findRotation(1, -2, -3, -4, 5, -6);

	arma::mat test1 = trans(rotated1) * rotated1;
	arma::mat test2 = trans(rotated2) * rotated2;
	arma::mat test3 = trans(rotated3) * rotated3;	
	arma::mat test4 = trans(rotated4) * rotated4;
	for(int k)
	return 0;
}

/** \brief Test for huckel theory initial guess.*/
int run_huckel_test(){ 
	arma::mat SMatrix;
    SMatrix.load("s.txt", arma::raw_ascii);

    std::vector<Semi::xyz> xyzData;
    std::ifstream fin;
    std::string line;
    int i = 0;
    fin.open("benzene.txt");
    while (std::getline(fin, line)) {
        std::stringstream linestream(line);
        double elem, x, y, z;
        linestream >> elem >> x >> y >> z;
        xyzData.push_back(Semi::xyz());
        xyzData[i].atom = elem;
        xyzData[i].x = x;
        xyzData[i].y = y;
        xyzData[i].z = z;
        i++;
    }
    fin.close();

    Semi::huckel(SMatrix, 1, 0.2, xyzData);
    return 0;
}

/** Main method to run tests.
 */
int main() {
	double result = run_compilation_test() | run_sto_test() | run_rotation_test() | 0;
	return result;
}
