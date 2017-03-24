#define ARMA_DONT_USE_WRAPPER
#include <iostream>
#include <vector>
#include <armadillo>
#include <semi/Atom.h>
#include <semi/Molecule.h>
#include <semi/Basis.h>
#include <semi/BasisSet.h>
#include "IntegralEvaluator.h"
using namespace arma;

namespace Semi {
bool isReversed(int *a, int *b) {
	if ((b[0] == 0 && b[1] == 0 && a[0] == 1 && a[1] == 0) ||
	        (b[0] == 1 && b[1] == 0 && a[0] == 1 && a[1] == 0) ||
	        (b[0] == 0 && b[1] == 0 && a[0] == 2 && a[1] == 0) ||
	        (b[0] == 1 && b[1] == 0 && a[0] == 2 && a[1] == 0) ||
	        (b[0] == 2 && b[1] == 0 && a[0] == 2 && a[1] == 0) ||
	        (b[0] == 0 && b[1] == 0 && a[0] == 2 && a[1] == 1 && a[2] == b[2]) ||
	        (b[0] == 1 && b[1] == 0 && a[0] == 2 && a[1] == 1 && a[2] == b[2]) ||
	        (b[0] == 2 && b[1] == 0 && a[0] == 2 && a[1] == 1 && a[2] == b[2]) ||
	        (b[0] == 1 && b[1] == 0 && a[0] == 1 && a[1] == 1 && a[2] == b[2]) ||
	        (b[0] == 1 && b[1] == 1 && a[0] == 2 && a[1] == 0 && a[2] == b[2]) ||
	        (b[0] == 1 && b[1] == 1 && a[0] == 2 && a[1] == 1 && a[2] == b[2]) ||
	        (b[0] == 2 && b[1] == 1 && a[0] == 2 && a[1] == 1 && a[2] == b[2]) ||
	        (b[0] == 1 && b[1] == 1 && a[0] == 2 && a[1] == 1 && a[2] == b[2]) ||
	        (b[0] == 2 && b[1] == 1 && a[0] == 2 && a[1] == 1 && a[2] == b[2])) {
		return true;
	}
	return false;
}

double CalculateOverlapFull(double tau, double rho, double kappa, double rho_alpha, double rho_beta, int *a, int *b) {
	double overlap = 0.0;

	if (a[0] == 0 && a[1] == 0 && b[0] == 1 && b[1] == 0) { //0s1s
		overlap = (pow((1.0 - pow(tau, 2)), 0.5) / (pow(2, 0.5) * tau * rho))
		          * (-(1.0 - kappa) * exp(-rho_alpha) + ((1.0 - kappa) + rho_beta) * exp(-rho_beta));
		return overlap;
	}
	if (a[0] == 1 && a[1] == 0 && b[0] == 1 && b[1] == 0) { //1s1s

		overlap = (pow((1.0 - pow(tau, 2)), 0.5) / (tau * rho))
		          * (-(1.0 - kappa) * (2.0 * (1.0 + kappa) + rho_alpha) * exp(-rho_alpha)
		             + (1.0 + kappa) * (2.0 * (1.0 - kappa) + rho_beta) * exp(-rho_beta));
		return overlap;
	}
	if (a[0] == 0 && a[1] == 0 && b[0] == 2 && b[1] == 0) { //0s2s
		overlap = (pow((1.0 - pow(tau, 2)), 0.5) / (pow(6, 0.5) * tau * rho))
		          * (-(1.0 - kappa) * (1.0 - 2.0 * kappa) * exp(-rho_alpha)
		             + ((1.0 - kappa) * (1.0 - 2.0 * kappa) + 2.0 * (1.0 - kappa) * rho_beta + pow(rho_beta, 2)) * exp(-rho_beta));
		return overlap;
	}
	if (a[0] == 1 && a[1] == 0 && b[0] == 2 && b[1] == 0) { //1s2s
		overlap = (pow((1.0 - pow(tau, 2)), 0.5) / (pow(3, 0.5) * tau * rho))
		          * (-(1.0 - kappa) * (2.0 * (1.0 + kappa) * (2.0 - 3.0 * kappa) + (1.0 - 2.0 * kappa) * exp(-rho_alpha))
		             + (1.0 + kappa) * (2.0 * (1.0 - kappa) * (2.0 - 3.0 * kappa) + 4.0 * (1.0 - kappa) * rho_beta + pow(rho_beta, 2)) * exp(-rho_beta));
		return overlap;
	}
	if (a[0] == 2 && a[1] == 0 && b[0] == 2 && b[1] == 0) { //2s2s
		overlap = (pow((1.0 - pow(tau, 2)), 0.5) / (3.0 * tau * rho))
		          * (-(1.0 - kappa) * (2.0 * (1.0 + kappa) * (7 - 12.0 * pow(kappa, 2)) + 4.0 * (1.0 + kappa) * (2.0 - 3.0 * kappa) * rho_alpha + (1.0 - 2.0 * kappa) * pow(rho_alpha, 2) *  exp(-rho_alpha))
		             + (1.0 + kappa) * (2.0 * (1.0 - kappa) * (7 - 12.0 * pow(kappa, 2)) + 4.0 * (1.0 - kappa) * (2.0 + 3.0 * kappa) * rho_beta + (1.0 + 2.0 * kappa) * pow(rho_beta, 2) *  exp(-rho_beta)));
		return overlap;
	}
	if (a[0] == 0 && a[1] == 0 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //0s2pz
		overlap = pow((1.0 + tau) / (1.0 - tau), 0.5) * (1.0 / (pow(2, 0.5) * tau * pow(rho, 2)))
		          * (-2.0 * pow(1.0 - kappa, 2) * (1.0 + rho_alpha) * exp(-rho_alpha)
		             + (2.0 * pow(1.0 - kappa, 2) * (1.0 + rho_beta) + 2.0 * (1.0 - kappa) * pow(rho_beta, 2) + pow(rho_beta, 3)) * exp(-rho_beta));
		return overlap;
	}
	if (a[0] == 1 && a[1] == 0 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //1s2pz
		overlap = pow((1.0 + tau) / (1.0 - tau), 0.5) * (1.0 / (tau * pow(rho, 2)))
		          * (-pow(1.0 - kappa, 2) * (6.0 * (1.0 + kappa) * (1.0 + rho_alpha) + 2.0 * pow(rho_alpha, 2)) * exp(-rho_alpha)
		             + pow(1.0 + kappa, 2) * (6.0 * (1.0 - kappa) * (1.0 + rho_beta) + 4.0 * (1.0 - kappa) * pow(rho_beta, 2) + pow(rho_beta, 3)) * exp(-rho_beta));
		return overlap;
	}
	if (a[0] == 2 && a[1] == 0 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //2s2pz
		overlap = pow((1.0 + tau) / (1.0 - tau), 0.5) * (1.0 / (pow(3, 0.5) * tau * pow(rho, 2)))
		          * (-pow(1.0 - kappa, 2) * (6.0 * (1.0 + kappa) * (3.0 + 4.0 * kappa) * (1.0 + rho_alpha) + 2.0 * (5.0 + 6.0 * kappa) * pow(rho_alpha, 2) + 2.0 * pow(rho_alpha, 3)) * exp(-rho_alpha)
		             + (1.0 + kappa) * (6.0 * pow(1.0 - kappa, 2) * (3.0 + 4.0 * kappa) * (1.0 + rho_beta) + 4.0 * (1.0 - kappa) * (2.0 + 3.0 * kappa) * pow(rho_beta, 2) + (1.0 + 2.0 * kappa) * pow(rho_alpha, 3)) * exp(-rho_beta));
		return overlap;
	}
	if (a[0] == 1 && a[1] == 0 && b[0] == 1 && b[1] == 1 && a[2] == b[2]) { //1pz1s
		overlap = pow((1.0 - tau) / (1.0 + tau), 0.5) * (pow(3, 0.5) / (tau * pow(rho, 2)))
		          * (-(1.0 - kappa) * (2.0 * (1.0 + kappa) * (1.0 + rho_alpha) + pow(rho_alpha, 2)) * exp(-rho_alpha)
		             + (1.0 + kappa) * (2.0 * (1.0 - kappa) * (1.0 + rho_beta) + pow(rho_beta, 2)) * exp(-rho_beta));
		return overlap;
	}
	if (a[0] == 1 && a[1] == 1 && b[0] == 2 && b[1] == 0 && a[2] == b[2]) { //1pz2s
		overlap = pow((1.0 - tau) / (1.0 + tau), 0.5) * (1.0 / (tau * pow(rho, 2)))
		          * (-(1.0 - kappa) * (2.0 * (1.0 + kappa) * (2.0 - 3.0 * kappa) * (1.0 + rho_alpha) + (1.0 - 2.0 * kappa) * pow(rho_alpha, 2)) * exp(-rho_alpha)
		             + (1.0 + kappa) * (2.0 * (1.0 - kappa) * (2.0 - 3.0 * kappa) * (1.0 + rho_beta) + (3.0 - 4.0 * kappa) * pow(rho_beta, 2) + pow(rho_beta, 3)) * exp(-rho_beta));
		return overlap;
	}
	if (a[0] == 1 && a[1] == 1 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //1pz2pz
		overlap = pow(3, 0.5) / (pow(1.0 - pow(tau, 2), 0.5) * tau * pow(rho, 3))
		          * (-pow(1.0 - kappa, 2) * (12.0 * (1.0 + kappa) * (1.0 + rho_alpha + 0.5 * pow(rho_alpha, 2)) + 2.0 * pow(rho_alpha, 3)) * exp(-rho_alpha)
		             + (1.0 + kappa) * (12.0 * pow(1.0 - kappa, 2) * (1.0 + rho_beta + 0.5 * pow(rho_beta, 2)) + (3.0 - 4.0 * kappa) * pow(rho_beta, 3) + pow(rho_beta, 4)) * exp(-rho_beta));
		return overlap;
	}
	if (a[0] == 2 && a[1] == 1 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //2pz2pz
		overlap = 1.0 / (pow(1.0 - pow(tau, 2), 0.5) * tau * pow(rho, 3))
		          * (-pow(1.0 - kappa, 2) * (48.0 * pow(1.0 + kappa, 2) * (1.0 + rho_alpha + 0.5 * pow(rho_alpha, 2)) + 2.0 * (5.0 + 6.0 * kappa) * pow(rho_alpha, 3) + 2.0 * pow(rho_alpha, 4)) * exp(-rho_alpha)
		             + pow(1.0 + kappa, 2) * (48.0 * pow(1.0 - kappa, 2) * (1.0 + rho_beta + 0.5 * pow(rho_beta, 2)) + 2.0 * (5.0 - 6.0 * kappa) * pow(rho_beta, 3) + 2.0 * pow(rho_beta, 4)) * exp(-rho_beta));
		return overlap;
	}
	if (a[0] == 1 && a[1] == 1 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //1pz2pz
		overlap = pow(3, 0.5) / (pow(1.0 - pow(tau, 2), 0.5) * tau * pow(rho, 3))
		          * (-pow(1.0 - kappa, 2) * (6.0 * (1.0 + kappa) * (1.0 + rho_alpha) + 2.0 * pow(rho_alpha, 2)) * exp(-rho_alpha)
		             + (1.0 + kappa) * (6.0 * pow(1.0 - kappa, 2) * (1.0 + rho_beta) + 4.0 * (1.0 - kappa) * pow(rho_beta, 2) + pow(rho_beta, 3)) * exp(-rho_beta));
		return overlap;
	}
	if (a[0] == 2 && a[1] == 1 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //2pz2pz
		overlap = 1.0 / (pow(1.0 - pow(tau, 2), 0.5) * tau * pow(rho, 3))
		          * (-pow(1.0 - kappa, 2) * (24.0 * pow(1.0 + kappa, 2) * (1.0 + rho_alpha) + 12.0 * (1.0 + kappa) * pow(rho_alpha, 2) + pow(rho_alpha, 3)) * exp(-rho_alpha)
		             + pow(1.0 + kappa, 2) * (24.0 * pow(1.0 - kappa, 2) * (1.0 + rho_beta) + 12.0 * (1.0 - kappa) * pow(rho_beta, 2) + pow(rho_beta, 3)) * exp(-rho_beta));
		return overlap;
	}
	if (isReversed(a, b)) {
		return CalculateOverlapFull(-tau, rho, -kappa, rho_alpha, rho_beta, a, b);
	}
}

double CalculateOverlapSameZeta(double tau, double rho, double kappa, double rho_alpha, double rho_beta, int *a, int *b) { //Tau = 0
	double overlap = 0.0;
	if (a[0] == 0 && a[1] == 0 && b[0] == 1 && b[1] == 0) { //0s1s
		overlap = (1.0 / pow(2, 0.5)) * (1.0 + rho) * exp(-rho);
		return overlap;
	}
	if (a[0] == 1 && a[1] == 0 && b[0] == 1 && b[1] == 0) { //1s1s
		overlap = (1.0 + rho + 1.0 / 3.0 * pow(rho, 2)) * exp(-rho);
		return overlap;
	}
	if (a[0] == 0 && a[1] == 0 && b[0] == 2 && b[1] == 0) { //0s2s
		overlap = pow(6, -0.5) * (1.0 + rho + 2.0 / 3.0 * pow(rho, 2)) * exp(-rho);
		return overlap;
	}
	if (a[0] == 1 && a[1] == 0 && b[0] == 2 && b[1] == 0) { //1s2s
		overlap = pow(3, 0.5) / 2.0 * (1.0 + rho + 4 / 9 * pow(rho, 2) + 1.0 / 9 * pow(rho, 3)) * exp(-rho);
		return overlap;
	}
	if (a[0] == 2 && a[1] == 0 && b[0] == 2 && b[1] == 0) { //2s2s
		overlap = (1.0 + rho + 4 / 9 * pow(rho, 2) + 1.0 / 9 * pow(rho, 3) + 1.0 / 45 * pow(rho, 4)) * exp(-rho);
		return overlap;
	}
	if (a[0] == 0 && a[1] == 0 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //0s2pz
		overlap = (1.0 / 3.0 * pow(2, 0.5) * rho) * (1.0 + rho) * exp(-rho);
		return overlap;
	}
	if (a[0] == 1 && a[1] == 0 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //1s2pz
		overlap = 1.0 / 2.0 * rho * (1.0 + rho + 1.0 / 3.0 * pow(rho, 2)) * exp(-rho);
		return overlap;
	}
	if (a[0] == 2 && a[1] == 0 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //2s2pz
		overlap = 1.0 / (2.0 * pow(3, 0.5)) * rho * (1.0 + rho + 7 / 15 * pow(rho, 2) + 2.0 / 15 * pow(rho, 3)) * exp(-rho);
		return overlap;
	}
	if (a[0] == 1 && a[1] == 0 && b[0] == 1 && b[1] == 1 && a[2] == b[2]) { //1pz1s
		overlap = 1.0 / pow(3, 0.5) * rho * (1.0 + rho) * exp(-rho);
		return overlap;
	}
	if (a[0] == 1 && a[1] == 1 && b[0] == 2 && b[1] == 0 && a[2] == b[2]) { //1pz2s
		overlap = 1.0 / 6.0 * (1.0 + rho + pow(rho, 2)) * exp(-rho);
		return overlap;
	}
	if (a[0] == 1 && a[1] == 1 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //1pz2pz
		overlap = 1.0 / 2.0 * pow(3, 0.5) * (-1.0 - rho + 1.0 / 3.0 * pow(rho, 3)) * exp(-rho);
		return overlap;
	}
	if (a[0] == 2 && a[1] == 1 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //2pz2pz
		overlap = (-1.0 - rho - 1.0 / 5 * pow(rho, 2) + 2.0 / 15 * pow(rho, 3) + 1.0 / 15 * pow(rho, 4)) * exp(-rho);
		return overlap;
	}
	if (a[0] == 1 && a[1] == 1 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //1pz2pz
		overlap = 1.0 / 2.0 *  pow(3, 0.5) * (1.0 + rho + 1.0 / 3.0 * pow(rho, 2)) * exp(-rho);
		return overlap;
	}
	if (a[0] == 2 && a[1] == 1 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //2pz2pz
		overlap = (1.0 + rho + 2.0 / 5 * pow(rho, 2) + 1.0 / 15 * pow(rho, 3)) * exp(-rho);
		return overlap;
	}
	if (isReversed(a, b) && overlap != 0) {
		return CalculateOverlapSameZeta(-tau, rho, -kappa, rho_alpha, rho_beta, a, b);
	}
}

double CalculateOverlapSamePosition(double tau, double rho, double kappa, double rho_alpha, double rho_beta, int *a, int *b) { //Rho = 0
	double overlap = 0.0;
	if (a[0] == 0 && a[1] == 0 && b[0] == 1 && b[1] == 0) { //0s1s
		overlap = 1.0 / pow(2, 0.5) * pow(1.0 + rho, 0.5) * pow(1.0 - rho, 1.5);
		return overlap;
	}
	if (a[0] == 1 && a[1] == 0 && b[0] == 1 && b[1] == 0) { //1s1s
		overlap = pow(1.0 + rho, 1.5) * pow(1.0 - rho, 1.5);
		return overlap;
	}
	if (a[0] == 0 && a[1] == 0 && b[0] == 2 && b[1] == 0) { //0s2s
		overlap = pow(6, -0.5) * pow(1.0 + rho, 0.5) * pow(1.0 - rho, 2.5);
		return overlap;
	}
	if (a[0] == 1 && a[1] == 0 && b[0] == 2 && b[1] == 0) { //1s2s
		overlap = 1.0 / 2.0 * pow(3, 0.5) * pow(1.0 + rho, 1.5) * pow(1.0 - rho, 2.5);
		return overlap;
	}
	if (a[0] == 2 && a[1] == 0 && b[0] == 2 && b[1] == 0) { //2s2s
		overlap = pow(1.0 + rho, 2.5) * pow(1.0 - rho, 2.5);
		return overlap;
	}
	if (a[0] == 0 && a[1] == 0 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //0s2pz
		overlap = 0;
		return overlap;
	}
	if (a[0] == 1 && a[1] == 0 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //1s2pz
		overlap = 0;
		return overlap;
	}
	if (a[0] == 2 && a[1] == 0 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //2s2pz
		overlap = 0;
		return overlap;
	}
	if (a[0] == 1 && a[1] == 0 && b[0] == 1 && b[1] == 1 && a[2] == b[2]) { //1pz1s
		overlap = 0;
		return overlap;
	}
	if (a[0] == 1 && a[1] == 1 && b[0] == 2 && b[1] == 0 && a[2] == b[2]) { //1pz2s
		overlap = 0;
		return overlap;
	}
	if (a[0] == 1 && a[1] == 1 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //1pz2pz
		overlap = 1.0 / 2.0 * pow(3, 0.5) * pow(1.0 + rho, 1.5) * pow(1.0 - rho, 2.5);
		return overlap;
	}
	if (a[0] == 2 && a[1] == 1 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //2pz2pz
		return overlap;
		overlap = 1.0 / 2.0 * pow(3, 0.5) * pow(1.0 + rho, 1.5) * pow(1.0 - rho, 2.5);
	}
	return overlap;
	if (a[0] == 1 && a[1] == 1 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //1pz2pz
		overlap = pow(1.0 + rho, 2.5) * pow(1.0 - rho, 2.5);
	}
	return overlap;
	if (a[0] == 2 && a[1] == 1 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //2pz2pz
		overlap = pow(1.0 + rho, 2.5) * pow(1.0 - rho, 2.5);
		return overlap;
	}
	if (isReversed(a, b)) {
		return CalculateOverlapSamePosition(-tau, rho, -kappa, rho_alpha, rho_beta, a, b);
	}
}

double CalculateOverlap(double tau, double rho, double kappa, double rho_alpha, double rho_beta, int *a, int *b) {
	double overlap = 0.0;
	if (std::abs(tau) < tolerance) {
		overlap = CalculateOverlapSameZeta(tau, rho, kappa, rho_alpha, rho_beta, a, b);
	}
	else if (std::abs(rho) < tolerance) {
		overlap = CalculateOverlapSamePosition(tau, rho, kappa, rho_alpha, rho_beta, a, b);
	}
	else {
		overlap = CalculateOverlapFull(tau, rho, kappa, rho_alpha, rho_beta, a, b);
	}
	return overlap;
}

arma::mat findRotation(double x1, double y1, double z1, double x2, double y2, double z2) {
	double x = x2 - x1;
	double y = y2 - y1;
	double z = z2 - z1;

	//xz plane rotation
	mat xzRotation(3, 3);
	xzRotation(0, 0) = x / pow(pow(x, 2) + pow(y, 2), 0.5);
	xzRotation(0, 1) = y / pow(pow(x, 2) + pow(y, 2), 0.5);
	xzRotation(0, 2) = 0;
	xzRotation(1, 0) = -y / pow(pow(x, 2) + pow(y, 2), 0.5);
	xzRotation(1, 1) = x / pow(pow(x, 2) + pow(y, 2), 0.5);
	xzRotation(1, 2) = 0;
	xzRotation(2, 0) = 0;
	xzRotation(2, 1) = 0;
	xzRotation(2, 2) = 1;

	//rotation to z axis
	mat zRotation(3, 3);
	zRotation(0, 0) = z / pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 0.5);
	zRotation(0, 1) = 0;
	zRotation(0, 2) = (-pow(pow(x, 2) + pow(y, 2), 0.5)) / (pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 0.5));
	zRotation(1, 0) = 0;
	zRotation(1, 1) = 1;
	zRotation(1, 2) = 0;
	zRotation(2, 0) = (pow(pow(x, 2) + pow(y, 2), 0.5)) / (pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 0.5));
	zRotation(2, 1) = 0;
	zRotation(2, 2) = z / pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 0.5);

	mat rotationMatrix(3, 3);
	rotationMatrix = zRotation * xzRotation;
	return rotationMatrix;
}

int factorial(int n) {
	return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

} //namespace Semi