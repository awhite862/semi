#ifndef HUCKELBASIS_H
#define HUCKELBASIS_H
#include <armadillo>

namespace Semi {
/** \brief Class that solved does basis change from Huckel Basis to full Basis. **/

void constructCBasis(arma::mat c_full, arma::mat Smatrix, arma::mat &sol);

} // namespace Semi

#endif // HUCKELBASIS_H
