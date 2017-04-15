/** \brief A test class for semi.*/
#include <armadillo>
#include <cstdlib>
#include <fstream>
#include <semi/Integral/IntegralEvaluator.h>
#include <semi/Integral/Cndo.h>
#include <semi/Huckel/huckel.h>
#include <semi/semi_utils.h>
#include <map>

using namespace Semi;
/** \brief Basic compilation test for semi classes.*/
int run_compilation_test() {
/*    Atom* a = new Atom(1, 1, 1, 1, 1);
    Molecule* m = new Molecule(std::vector<Semi::Atom> (4, *a));

    arma::colvec r;
    r << 1 << 1 << 1;
    QNumber* q = new QNumber(1, 1, 0);
    STOBasis* sto = new STOBasis(1, 1, 1, 1, 1, 1, 1, 1);
    GTOBasis* gto = new GTOBasis(1, 1, 1, 1, r);
    //CGTOBasis* cgto = new CGTOBasis(0, 0, 0, r, 1, *q);*/

    return 0;
}

/** \brief Basic Overlap Integral test.
 *  Overlap of 2 1s orbitals, centered at (0, 0, 0), (1, 0, 0)
 *  1st test zeta values 1, 1, evaluates to 0.214596340683341354 un-normalized
 *  2nd test zeta values 0.2, 1, evaluates to 0.121087705563123644 un-normalized
 */
int run_sto_test() {
    ///Unnormalized values
    double actualFirst = 0.214596340683341354;
    double actualSecond = 1.121087705563123644;

    ///paramaters for the first overlap integral
    double r = 1, z1 = 1, z2 = 1;
    double tau = (z1 - z2) / (z1 + z2);
    double rho = 0.5 * (z1 + z2) * r;
    double kappa = 0.5 * (tau + 1 / tau);
    double rho_alpha = z1 * r, rho_beta = z2 * r;
    int a[3] = {1, 0, 0};
    int b[3] = {1, 0, 0};
    double first = calculateOverlapSTO(tau, rho, kappa, rho_alpha, rho_beta, a, b);
    double normFirst = pow(2 * z1, 1 + 0.5) * pow(Semi::factorial(2 * 1), -0.5) *  pow(2 * z2, 1 + 0.5) * pow(Semi::factorial(2 * 1), -0.5) * actualFirst;

    ///parameters for second overlap integral
    z1 = 0.2; z2 = 1;
    tau = (z1 - z2) / (z1 + z2);
    rho = 0.5 * (z1 + z2) * r;
    kappa = 0.5 * (tau + 1.0 / tau);
    rho_alpha = z1 * r; rho_beta = z2 * r;
    double second = Semi::calculateOverlapSTO(tau, rho, kappa, rho_alpha, rho_beta, a, b);
    double normSecond = pow(2 * z1, 1 + 0.5) * pow(Semi::factorial(2 * 1), -0.5) *  pow(2 * z2, 1 + 0.5) * pow(Semi::factorial(2 * 1), -0.5) * actualSecond;

    return !(abs(normFirst-first) < tolerance || abs(normSecond-second) < tolerance);
}

/** \brief Test for rotation matrix in overlap integrals.*/
int run_rotation_test() {
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

    arma::vec vec5(3);
    vec4(0) = 0;
    vec4(1) = 0;
    vec4(2) = 5;
    arma::mat rotated5 = Semi::findRotation(0, 0, 5, 0, 0, 10);

    if (norm(test1 - arma::eye(3, 3)) > tolerance | norm(test2 - arma::eye(3, 3)) > tolerance | norm(test3 - arma::eye(3, 3)) > tolerance | norm(test4 - arma::eye(3, 3)) > tolerance | norm(rotated5 - arma::eye(3, 3)) > tolerance) {
        return 1;
    }

    return 0;
}

/** \brief Test for huckel theory initial guess.*/
int run_huckel_test() {
    arma::mat SMatrix;
    SMatrix.load("s.txt", arma::raw_ascii);

    Molecule m;
    std::ifstream fin;
    std::string line;
    int i = 0;
    fin.open("benzene.txt");
    while (std::getline(fin, line)) {
        std::stringstream linestream(line);
        double elem, x, y, z;
        linestream >> elem >> x >> y >> z;
        Atom* a = new Atom(x,y,z,elem,i);
        m.myMolecule.push_back(*a);
        i++;
    }
    fin.close();

    Semi::calculateHuckel(SMatrix, 1, 0.2, m, "c_v");
    return 0;
}
// 0 0 l

// 1.0 0.01 0.001 abc
// -1.0 -0.01 -0.001 abc

// 10.0 10.0

// 4.007837386698362180096159730696235475449195275225351306619521440472961255744995533476335816628324
// 4.007837386698362180096159730696235475449195275225351306619521440472961255744995533476335816628324

// 2.056994294456803574829111780821946482718902295308792506777564143425656498691520607978489E-9
/** \brief Test CGTO overlap.*/
int run_cgto_test() {
    // arma::colvec ar;
    // ar << 1 << 0.01 << 0.001;
    // QNumber* aq = new QNumber(1, 0, 0);
    // CGTOBasis a = CGTOBasis(*aq, 0, 0, 0, ar, 1);

    // arma::colvec br;
    // br << -1 << -0.01 << -0.001;
    // QNumber* bq = new QNumber(1, 0, 0);
    // CGTOBasis b = CGTOBasis(*bq, 0, 0, 0, br, 1);

    // double ans = calculateOverlapCGTO(a, b);
    arma::colvec ar;
    ar << 1 << 0.01 << 0.001;
    QNumber* aq = new QNumber(1, 0, 0);
    GTOBasis a = GTOBasis(*aq, 0, 0, 0, 10, ar);

    arma::colvec br;
    br << -1 << -0.01 << -0.001;
    QNumber* bq = new QNumber(1, 0, 0);
    GTOBasis b = GTOBasis(*bq, 0, 0, 0, 10, br);

    double ans = calculateOverlapGTO(a, b);
    
    arma::colvec ar2;
    ar2 << 1 << 0.01 << 0.001;
    QNumber* aq2 = new QNumber(1, 1, 0);
    GTOBasis a2 = GTOBasis(*aq2, 1, 0, 0, 10, ar2);

    arma::colvec br2;
    br2 << -1 << -0.01 << -0.001;
    QNumber* bq2 = new QNumber(1, 1, 0);
    GTOBasis b2 = GTOBasis(*bq2, 1, 0, 0, 10, br2);

    double ans2 = calculateOverlapGTO(a2, b2);
    
    arma::colvec ar3;
    ar3 << 10 << 1 << 0.1;
    QNumber* aq3 = new QNumber(1, 1, 0);
    GTOBasis a3 = GTOBasis(*aq3, 1, 0, 0, 0.01, ar3);

    arma::colvec br3;
    br3 << -10 << -1 << -0.1;
    QNumber* bq3 = new QNumber(1, 1, 0);
    GTOBasis b3 = GTOBasis(*bq3, 1, 0, 0, 0.05, br3);

    double ans3  = calculateOverlapGTO(a3, b3);
    std::cout << ans << std::endl;
    std::cout << ans2 << std::endl;
    std::cout << ans3 << std::endl;
    return 0;
}

int run_cgto_matrix_test(){
    arma::colvec ar;
    ar << 1 << 0.01 << 0.001;
    QNumber* aq = new QNumber(1, 0, 0);
    GTOBasis a = GTOBasis(*aq, 0, 0, 0, 3.4252509, ar);

    arma::colvec br;
    br << 1 << 0.01 << 0.001;
    QNumber* bq = new QNumber(1, 0, 0);
    GTOBasis b = GTOBasis(*bq, 0, 0, 0, 0.6239137, br);        

    arma::colvec cr;
    cr << 1 << 0.01 << 0.001;
    QNumber* cq = new QNumber(1, 0, 0);
    GTOBasis c = GTOBasis(*cq, 0, 0, 0, 0.1688554, cr);

    arma::colvec dr;
    dr << 2 << 0.01 << 0.001;
    QNumber* dq = new QNumber(1, 0, 0);
    GTOBasis d = GTOBasis(*dq, 0, 0, 0, 3.4252509, dr);

    arma::colvec er;
    er << 2 << 0.01 << 0.001;
    QNumber* eq = new QNumber(1, 0, 0);
    GTOBasis e = GTOBasis(*eq, 0, 0, 0, 0.6239137, er);


    arma::colvec fr;
    fr << 2 << 0.01 << 0.001;
    QNumber* fq = new QNumber(1, 0, 0);
    GTOBasis f = GTOBasis(*fq, 0, 0, 0, 0.1688554, fr);

    CGTOBasis abc = CGTOBasis(*aq, 0,0,0, ar, 1);

    CGTOBasis def = CGTOBasis(*dq, 0,0,0, dr, 1);
    double overlap = calculateOverlapCGTO(abc, def);
    double overlap_sum = 0
        + calculateOverlapGTO(a,d) * 0.15432897 * 0.15432897
        + calculateOverlapGTO(a,e) * 0.15432897 * 0.53532814
        + calculateOverlapGTO(a,f) * 0.15432897 * 0.44463454
        + calculateOverlapGTO(b,d) * 0.53532814 * 0.15432897
        + calculateOverlapGTO(b,e) * 0.53532814 * 0.53532814
        + calculateOverlapGTO(b,f) * 0.53532814 * 0.44463454
        + calculateOverlapGTO(c,d) * 0.44463454 * 0.15432897
        + calculateOverlapGTO(c,e) * 0.44463454 * 0.53532814
        + calculateOverlapGTO(c,f) * 0.44463454 * 0.44463454;



        std::cout << "sum" << std::endl;
        std::cout << calculateOverlapGTO(a,d) * 0.15432897 * 0.15432897 << std::endl;
        std::cout << calculateOverlapGTO(a,e) * 0.15432897 * 0.53532814 << std::endl;
        std::cout << calculateOverlapGTO(a,f) * 0.15432897 * 0.44463454 << std::endl;
        std::cout << calculateOverlapGTO(b,d) * 0.53532814 * 0.15432897 << std::endl;
        std::cout << calculateOverlapGTO(b,e) * 0.53532814 * 0.53532814 << std::endl;
        std::cout << calculateOverlapGTO(b,f) * 0.53532814 * 0.44463454 << std::endl;
        std::cout << calculateOverlapGTO(c,d) * 0.44463454 * 0.15432897 << std::endl;
        std::cout << calculateOverlapGTO(c,e) * 0.44463454 * 0.53532814 << std::endl;
        std::cout << calculateOverlapGTO(c,f) * 0.44463454 * 0.44463454 << std::endl;







        std::cout << "norm: " << overlap << std::endl;
        std::cout << "sum: " << overlap_sum << std::endl;
    return 0;
}

/** Main method to run tests.
 */
int main() {
    double result = run_compilation_test() | run_sto_test() | run_rotation_test() | run_huckel_test() | run_cgto_matrix_test() | 0;
    return result;
}
