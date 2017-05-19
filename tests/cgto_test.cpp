/** \brief A test class for semi.*/
#include <armadillo>
#include <cstdlib>
#include <fstream>
#include <semi/Integral/IntegralEvaluator.h>
#include <semi/Integral/Cndo.h>
#include <semi/Huckel/huckel.h>
#include <semi/semi_utils.h>
#include <semi/CGTOBasis.h>
#include <semi/GTOBasis.h>
#include <semi/QNumber.h>
#include <map>

using namespace Semi;

/** \brief Test for evaluating simple GTO Overlaps.*/

int run_cgto_test() {
    arma::colvec ar;
    ar << 1 << 0.01 << 0.001;
    QNumber aq(1, 0, 0);
    GTOBasis a(aq, 0, 0, 0, 10, ar);

    arma::colvec br;
    br << -1 << -0.01 << -0.001;
    QNumber bq(1, 0, 0);
    GTOBasis b(bq, 0, 0, 0, 10, br);

    double ans1 = calculateOverlapGTO(a, b);

    arma::colvec ar2;
    ar2 << 1 << 0.01 << 0.001;
    QNumber aq2(1, 1, 0);
    GTOBasis a2(aq2, 1, 0, 0, 10, ar2);

    arma::colvec br2;
    br2 << -1 << -0.01 << -0.001;
    QNumber bq2(1, 1, 0);
    GTOBasis b2(bq2, 1, 0, 0, 10, br2);

    double ans2 = calculateOverlapGTO(a2, b2);

    arma::colvec ar3;
    ar3 << 10 << 1 << 0.1;
    QNumber aq3(1, 1, 0);
    GTOBasis a3(aq3, 1, 0, 0, 0.01, ar3);

    arma::colvec br3;
    br3 << -10 << -1 << -0.1;
    QNumber bq3(1, 1, 0);
    GTOBasis b3(bq3, 1, 0, 0, 0.05, br3);

    double ans3  = calculateOverlapGTO(a3, b3);
    double calc1 = 2.056994294456803574829111780821946482718902295308792506777564143425656498691520607978489 * pow(10, -9);
    double calc2 = -8.0222777483815339418335359452055912826037189517042907764325001593600603458514467053664595 * pow(10, -8);
    double calc3 = -0.093749187977773697056124931824321559965875946561514552449808814321843488815579167171142953324271;
    return abs(ans1 - calc1) > tolerance | abs(ans2 - calc2) > tolerance | abs(ans3 - calc3) > tolerance;
}

/** \brief Test for if CGTO from formula is the same as summing multiple GTOs.*/
int run_cgto_matrix_test() {
    arma::colvec ar;
    ar << 1 << 0.01 << 0.001;
    QNumber aq(1, 0, 0);
    GTOBasis a(aq, 0, 0, 0, 3.4252509, ar);

    arma::colvec br;
    br << 1 << 0.01 << 0.001;
    QNumber bq(1, 0, 0);
    GTOBasis b(bq, 0, 0, 0, 0.6239137, br);

    arma::colvec cr;
    cr << 1 << 0.01 << 0.001;
    QNumber cq(1, 0, 0);
    GTOBasis c(cq, 0, 0, 0, 0.1688554, cr);

    arma::colvec dr;
    dr << 2 << 0.01 << 0.001;
    QNumber dq(1, 0, 0);
    GTOBasis d(dq, 0, 0, 0, 3.4252509, dr);

    arma::colvec er;
    er << 2 << 0.01 << 0.001;
    QNumber eq(1, 0, 0);
    GTOBasis e(eq, 0, 0, 0, 0.6239137, er);

    arma::colvec fr;
    fr << 2 << 0.01 << 0.001;
    QNumber fq(1, 0, 0);
    GTOBasis f(fq, 0, 0, 0, 0.1688554, fr);

    CGTOBasis abc = CGTOBasis(aq, 0, 0, 0, ar, 1);

    CGTOBasis def = CGTOBasis(dq, 0, 0, 0, dr, 1);
    double overlap_calc = calculateOverlapCGTO(abc, def);
    double overlap_sum = 0
                         + calculateOverlapGTO(a, d) * 0.15432897 * 0.15432897
                         + calculateOverlapGTO(a, e) * 0.15432897 * 0.53532814
                         + calculateOverlapGTO(a, f) * 0.15432897 * 0.44463454
                         + calculateOverlapGTO(b, d) * 0.53532814 * 0.15432897
                         + calculateOverlapGTO(b, e) * 0.53532814 * 0.53532814
                         + calculateOverlapGTO(b, f) * 0.53532814 * 0.44463454
                         + calculateOverlapGTO(c, d) * 0.44463454 * 0.15432897
                         + calculateOverlapGTO(c, e) * 0.44463454 * 0.53532814
                         + calculateOverlapGTO(c, f) * 0.44463454 * 0.44463454;

    return abs(overlap_calc - overlap_sum) > tolerance;
}

/** Main method to run tests.
 */
int main() {
    double result = run_cgto_test() | run_cgto_matrix_test() | 0;
    return result;
}
