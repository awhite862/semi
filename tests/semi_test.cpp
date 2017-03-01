#include <fstream>
#include <io/io_input_section.h>
#include <cstdlib> // EXIT_FAILURE failed compilation
#include "../src/semi/Atom.h"
#include "../src/semi/Basis.h"
#include "../src/semi/Integral/IntegralEvaluator.h"
#include "../src/semi/Integral/Cndo.h"
#include <armadillo>

/** \test basic compilation. **/
int run_compilation_test() {
	Semi::Atom* a = new Semi::Atom(1, 1, 1, 1);
	Semi::Basis* b = new Semi::Basis(1, 1, 1, 1, 1, 1, 1);
	return 0;
}


//centers 0 0 0, 1 0 0
//exponents 1 1 (integral 1) 0.2 1 (integral 2)
//0.214596340683341354       1.121087705563123644
//1s functions


int run_sto_test() {
	double r = 1;
	double z1 = 1;
	double z2 = 1;
	double tau = (z1 - z2) / (z1 + z2);
	double rho = 0.5 * (z1 + z2) * r;
	double kappa = 0.5 * (rho + 1 / rho);
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
	double first = Semi::CalculateOverlap(tau, rho, kappa, rho_alpha, rho_beta, a, b);

	z1 = 0.2;
	z2 = 1;
	tau = (z1 - z2) / (z1 + z2);
	rho = 0.5 * (z1 + z2) * r;
	kappa = 0.5 * (rho + 1 / rho);
	rho_alpha = z1 * r;
	rho_beta = z2 * r;
	double second = Semi::CalculateOverlap(tau, rho, kappa, rho_alpha, rho_beta, a, b);

	double actualFirst = 0.214596340683341354;
	double actualSecond = 1.121087705563123644;

	first -= actualFirst;
	second -= actualSecond;

	return abs(second) < 0.00000000000001;
}

int main() {
	std::ostringstream oss;
    std::streambuf* p_cout_streambuf = std::cout.rdbuf();
    std::cout.rdbuf(oss.rdbuf());

	double result = run_compilation_test() | run_sto_test() | 0;

	std::cout.rdbuf(p_cout_streambuf);
	std::cout << oss.str();

	return result;
}
