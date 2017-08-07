#include "CGTOFunction.h"
#include "get_sto_3g.h"
#include "semi/semi_utils.h"

namespace Semi {

CGTOFunction::CGTOFunction(QNumber _nlm, double _a, double _b, double _c, arma::colvec _r, double _elem) {
    a = _a;
    b = _b;
    c = _c;
    r = _r;
    elem = _elem;
    nlm = _nlm;

    //loads sto_3g contraction coefficients 
    std::string name = getElement(_elem);
    alphaVec.resize(3);
    nVec.resize(3);
    get_sto_3g(name, nlm, nVec, alphaVec);

    for (int k = 0; k < 3; k++) {
        double n = pow(pow(M_PI / (2 * alphaVec[k]), 3.0 / 2.0)
                       * (doubleFactorial(2.0 * a - 1) * doubleFactorial(2.0 * b - 1)
                          * doubleFactorial(2.0 * c - 1)) / (pow(2, 2.0 * nlm.l) * pow(alphaVec[k], nlm.l)), -1.0 / 2.0);
        nVec[k] *= n;
    }
}

} // namespace Semi
