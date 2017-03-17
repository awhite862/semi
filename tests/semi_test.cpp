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
#include <armadillo>
using namespace Semi;
/** \brief Basic compilation test for semi classes.*/
int run_compilation_test() {
	Semi::Atom * a = new Semi::Atom(1, 1, 1, 1);
	Semi::Basis * b = new Semi::Basis(1, 1, 1, 1, 1, 1, 1);
	Semi::Molecule * Molecule = new Semi::Molecule(std::vector<Semi::Atom> (4, *a));
	Semi::BasisSet * BasisSet = new Semi::BasisSet(std::vector<Semi::Basis> (4, *b));
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
	CalculateOverlap(tau, rho, kappa, rho_alpha, rho_beta, a, b);
	double normFirst = pow(2 * z1, 1 + 0.5) * pow(Semi::IntegralEvaluator().factorial(2 * 1), -0.5) *  pow(2 * z2, 1 + 0.5) * pow(Semi::IntegralEvaluator().factorial(2 * 1), -0.5) * actualFirst;

	///parameters for second overlap integral
	z1 = 0.2;
	z2 = 1;
	tau = (z1 - z2) / (z1 + z2);
	rho = 0.5 * (z1 + z2) * r;
	kappa = 0.5 * (tau + 1.0 / tau);
	rho_alpha = z1 * r;
	rho_beta = z2 * r;
	double second = Semi::IntegralEvaluator().CalculateOverlap(tau, rho, kappa, rho_alpha, rho_beta, a, b);
	double normSecond = pow(2 * z1, 1 + 0.5) * pow(Semi::IntegralEvaluator().factorial(2 * 1), -0.5) *  pow(2 * z2, 1 + 0.5) * pow(Semi::IntegralEvaluator().factorial(2 * 1), -0.5) * actualSecond;
	
	normFirst -= first;
	normSecond -= second;
	return abs(normFirst) > tol || abs(normSecond) > tol;
}

/** \brief Test for rotation matrix in overlap integrals.*/
int run_rotation_test() {
	arma::vec x(3);
	x(1) = 5;
	x(2) = 5;
	x(3) = 5;
	// arma::mat rotation = Semi::findRotation(5, 5, 5, 0, 0, 0);
	// arma::vec rotated = rotation * x;
	// rotated.print();
	return 0;
}

/** Main method to run tests.
 */
int main() {
	double result = run_compilation_test() | run_sto_test() | run_rotation_test() | 0;
	return result;
}
