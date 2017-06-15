#include "HuckelBasis.h"
#include "semi/semi_utils.h"
namespace Semi {

void constructCBasis(arma::mat c_full, arma::mat Smatrix, arma::mat &sol) {
	arma::mat y = trans(c_full) * Smatrix * c_full;
	y.print();
	arma::mat x;
	invSqrt(y, x);
	x.print();
	arma::mat temp;
	invSqrt(y, temp);
	arma::mat c_basis = c_full * temp;
	arma::mat id = trans(c_basis) * Smatrix * c_basis;
	c_basis.print("c_basis");
	sol = c_basis;
}

} // namespace Semi

