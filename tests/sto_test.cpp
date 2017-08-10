/** \brief A test class for semi.*/
#include <armadillo>
#include <cstdlib>
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

/** \brief Basic Overlap Integral test.
 *  Overlap of 2 1s orbitals, centered at (0, 0, 0), (1, 0, 0)
 *  1st test zeta values 1, 1, evaluates to 0.214596340683341354 un-normalized
 *  2nd test zeta values 0.2, 1, evaluates to 0.121087705563123644 un-normalized
 */
int run_sto_test() {
    double actualFirst = 0.214596340683341354;
    double actualSecond = 1.121087705563123644;

    double r = 1, z1 = 1, z2 = 1;
    double tau = (z1 - z2) / (z1 + z2);
    double rho = 0.5 * (z1 + z2) * r;
    double kappa = 0.5 * (tau + 1 / tau);
    double rho_alpha = z1 * r, rho_beta = z2 * r;

    QNumber aQNumber(1, 0, 0);
    QNumber bQNumber(1, 0, 0);
    double first = Semi::calculateOverlapSTO(tau, rho, kappa, rho_alpha, rho_beta, aQNumber, bQNumber);
    double normFirst = pow(2 * z1, 1 + 0.5) * pow(Semi::factorial(2 * 1), -0.5) *  pow(2 * z2, 1 + 0.5) * pow(Semi::factorial(2 * 1), -0.5) * actualFirst;

    z1 = 0.2; z2 = 1;
    tau = (z1 - z2) / (z1 + z2);
    rho = 0.5 * (z1 + z2) * r;
    kappa = 0.5 * (tau + 1.0 / tau);
    rho_alpha = z1 * r; rho_beta = z2 * r;

    double second = Semi::calculateOverlapSTO(tau, rho, kappa, rho_alpha, rho_beta, aQNumber, bQNumber);
    double normSecond = pow(2 * z1, 1 + 0.5) * pow(Semi::factorial(2 * 1), -0.5) *  pow(2 * z2, 1 + 0.5) * pow(Semi::factorial(2 * 1), -0.5) * actualSecond;

    return !(abs(normFirst - first) < tolerance || abs(normSecond - second) < tolerance);
}

int run_sto_matrix_test() {
    BasisSet<STOFunction> a;
    double f = 1.1823090398;
    QNumber q1s(1, 0, 0);
    QNumber q2s(2, 0, 0);
    QNumber q2px(2, 1, 1);
    QNumber q2py(2, 1, -1);
    QNumber q2pz(2, 1, 0);
    a.myBasis.push_back(STOFunction(q2s, 6, 0, 0, 0, 0));
    a.myBasis.push_back(STOFunction(q2px, 6, 0, 0, 0, 0));
    a.myBasis.push_back(STOFunction(q2py, 6, 0, 0, 0, 0));
    a.myBasis.push_back(STOFunction(q2pz, 6, 0, 0, 0, 0));

    a.myBasis.push_back(STOFunction(q1s, 1, f, f, f, 1));
    a.myBasis.push_back(STOFunction(q1s, 1, f, -f, -f, 2));
    a.myBasis.push_back(STOFunction(q1s, 1, -f, -f, f, 3));
    a.myBasis.push_back(STOFunction(q1s, 1, -f, f, -f, 4));
    arma::mat sol;
    Semi::calculateOverlapMatrix(a, sol);

    arma::mat eigvec;
    arma::vec eigval;
    eig_sym(eigval, eigvec, sol);

    arma::mat solution(8, 1);
    solution(0, 0) = 0.1994;
    solution(1, 0) = 0.3321;
    solution(2, 0) = 0.3321;
    solution(3, 0) = 0.3321;
    solution(4, 0) = 1.4802;
    solution(5, 0) = 1.4802;
    solution(6, 0) = 1.4802;
    solution(7, 0) = 2.3637;

    return norm(solution - eigval) > 1;
}

/** \brief Test for rotation matrix in overlap integrals.*/
int run_rotation_test() {
    arma::vec vec1(3);
    vec1(0) = 5;
    vec1(1) = 5;
    vec1(2) = 5;
    arma::mat rotated1;
    Semi::findRotation(5, 5, 5, 7, 7, 7, rotated1);

    arma::vec vec2(3);
    vec2(0) = 5;
    vec2(1) = 5;
    vec2(2) = 5;
    arma::mat rotated2;
    Semi::findRotation(5, 5, 5, 3, 3, 3, rotated2);

    arma::vec vec3(3);
    vec3(0) = 1;
    vec3(1) = 1;
    vec3(2) = 1;
    arma::mat rotated3;
    Semi::findRotation(1, 1, 1, 3, 4, 5, rotated3);

    arma::vec vec4(3);
    vec4(0) = 1;
    vec4(1) = -2;
    vec4(2) = -3;
    arma::mat rotated4;
    Semi::findRotation(1, -2, -3, -4, 5, -6, rotated4);

    arma::mat test1 = trans(rotated1) * rotated1;
    arma::mat test2 = trans(rotated2) * rotated2;
    arma::mat test3 = trans(rotated3) * rotated3;
    arma::mat test4 = trans(rotated4) * rotated4;

    arma::vec vec5(3);
    vec4(0) = 0;
    vec4(1) = 0;
    vec4(2) = 5;
    arma::mat rotated5;
    Semi::findRotation(0, 0, 5, 0, 0, 10, rotated5);

    return (norm(test1 - arma::eye(3, 3)) > tolerance | norm(test2 - arma::eye(3, 3)) > tolerance | norm(test3 - arma::eye(3, 3)) > tolerance | norm(test4 - arma::eye(3, 3)) > tolerance | norm(rotated5 - arma::eye(3, 3)) > tolerance);
}

/** Main method to run tests.
 */
int main() {
    double result = run_rotation_test() | run_sto_matrix_test() | run_rotation_test() | 0;
    return result;
}