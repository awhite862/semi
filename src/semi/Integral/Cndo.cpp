#define ARMA_DONT_USE_WRAPPER
#include <iostream>
#include <vector>
#include <armadillo>
#include "Atom.h"
#include "Molecule.h"
#include "Basis.h"
#include "BasisSet.h"
using namespace semi;
using namespace arma;

//U_uu
double calculateCoreHamiltonian() {

}

//P_AA
double calculateTotalChargeDensity() {

}

//P_uv
double calculateChargeDensity() {

}

//P_uv
double calculateTotalChargeDensity() {

}

//Gamma_uv
double calculateElectronRepulsion() {

}

//V_AB
double calculateNucleurAttraction() {

}

//B_AB
double calculateBondingParameter() {

}

//S_uv
double calculateOverlapMatrix(BasisSet a, BasisSet b) {
	mat Smatrix(a.size(), b.size());
	for (uint k = 0; k < a.size(); k++) {
		for (uint l = 0; l < b.size(); l++) {
			r = distance()
			    zeta_average = 0.5 * (a[k].zeta + b[l].zeta);
			tau = (a[k].zeta - b[l].zeta) / (a[k].zeta + b[l].zeta);
			rho = 0.5 * (a[k].zeta + b[l].zeta) * r;
			kappa = 0.5 * (rho + 1 / rho);
			rho_alpha = a[k].zeta * r;
			rho_beta = b[l].zeta * r;
			std::string overlap = a[k] + b[l];
			if (overlap == "0s1s") {
				Smatrix(k, l) = (pow((1 - pow * (tau, 2)), 0.5) / (pow(2, 0.5) * tau * rho))
				                * (-(1 - kappa) * exp(-rho_alpha) + ((1 - kappa) + rho_beta) * exp(-rho_beta));
			}
			if (overlap == "1s1s") {
				Smatrix(k, l) = (pow((1 - pow * (tau, 2)), 0.5) / (tau * rho))
				                * (-(1 - kappa) * (2 * (1 + kappa) + rho_alpha) * exp(-rho_alpha)
				                   + (1 + kappa) * (2 * (1 - kappa) + rho_beta) * exp(-rho_beta));
			}
			if (overlap == "0s2s") {
				Smatrix(k, l) = (pow((1 - pow * (tau, 2)), 0.5) / (pow(6, 0.5) * tau * rho))
				                * (-(1 - kappa) * (1 - 2 * kappa) * exp(-rho_alpha)
				                   + ((1 - kappa) * (1 - 2 * kappa) + 2 * (1 - kappa) * rho_beta + pow(rho_beta, 2)) * exp(-rho_beta));
			}
			if (overlap == "1s2s") {
				Smatrix(k, l) = (pow((1 - pow * (tau, 2)), 0.5) / (pow(3, 0.5) * tau * rho))
				                * (-(1 - kappa) * (2 * (1 + kappa) * (2 - 3 * kappa) + (1 - 2 * kappa) * exp(-rho_alpha))
				                   + (1 + kappa) * (2 * (1 - kappa) * (2 - 3 * kappa) + 4 * (1 - kappa) * rho_beta + pow(rho_beta, 2)) * exp(-rho_beta));
			}
			if (overlap == "2s2s") {
				s
				Smatrix(k, l) = (pow((1 - pow * (tau, 2)), 0.5) / (3 * tau * rho))
				                * (-(1 - kappa) * (2 * (1 + kappa) * (7 - 12 * pow(kappa, 2)) + 4 * (1 + kappa) * (2 - 3 * kappa) * rho_alpha + (1 - 2 * kappa) * pow(rho_alpha, 2) *  exp(-rho_alpha))
				                   + (1 + kappa) * (2 * (1 - kappa) * (7 - 12 * pow(kappa, 2)) + 4 * (1 - kappa) * (2 + 3 * kappa) * rho_beta + (1 + 2 * kappa) * pow(rho_beta, 2) *  exp(-rho_beta)));
			}
			if (overlap == "0s2pz") {
				Smatrix(k, l) = pow((1 + tau) / (1 - tau), 0.5) * (1 / (pow(2, 0.5) * tau * pow(rho, 2)))
				                * (-2 * pow(1 - kappa, 2) * (1 + rho_alpha) * exp(-rho_alpha)
				                   + (2 * pow(1 - kappa, 2) * (1 + rho_beta) + 2 * (1 - kappa) * pow(rho_beta, 2) + pow(rho_beta, 3)) * exp(-rho_beta));
			}
			if (overlap == "1s2pz") {
				Smatrix(k, l) = pow((1 + tau) / (1 - tau), 0.5) * (1 / (tau * pow(rho, 2)))
				                * (-pow(1 - kappa, 2) * (6 * (1 + kappa) * (1 + rho_alpha) + 2 * pow(rho_alpha, 2)) * exp(-rho_alpha)
				                   + pow(1 + kappa, 2) * (6 * (1 - kappa) * (1 + rho_beta) + 4 * (1 - kappa) * pow(rho_beta, 2) + pow(rho_beta, 3)) * exp(-rho_beta));
			}
			if (overlap == "2s2pz") {
				Smatrix(k, l) = pow((1 + tau) / (1 - tau), 0.5) * (1 / (pow(3, 0.5) * tau * pow(rho, 2)))
				                * (-pow(1 - kappa, 2) * (6 * (1 + kappa) * (3 + 4 * kappa) * (1 + rho_alpha) + 2 * (5 + 6 * kappa) * pow(rho_alpha, 2) + 2 * pow(rho_alpha, 3)) * exp(-rho_alpha)
				                   + (1 + kappa) * (6 * pow(1 - kappa, 2) * (3 + 4 * kappa) * (1 + rho_beta) + 4 * (1 - kappa) * (2 + 3 * kappa) * pow(rho_beta, 2) + (1 + 2 * kappa) * pow(rho_alpha, 3)) * exp(-rho_beta));
			}
			if (overlap == "1pz1s") {
				Smatrix(k, l) = pow((1 - tau) / (1 + tau), 0.5) * (pow(3, 0.5) / (tau * pow(rho, 2)))
				                * (-(1 - kappa) * (2 * (1 + kappa) * (1 + rho_alpha) + pow(rho_alpha, 2)) * exp(-rho_alpha)
				                   + (1 + kappa) * (2 * (1 - kappa) * (1 + rho_beta) + pow(rho_beta, 2)) * exp(-rho_beta));
			}
			if (overlap == "1pz2s") {
				Smatrix(k, l) = pow((1 - tau) / (1 + tau), 0.5) * (1 / (tau * pow(rho, 2)))
				                * (-(1 - kappa) * (2 * (1 + kappa) * (2 - 3 * kappa) * (1 + rho_alpha) + (1 - 2 * kappa) * pow(rho_alpha, 2)) * exp(-rho_alpha)
				                   + (1 + kappa) * (2 * (1 - kappa) * (2 - 3 * kappa) * (1 + rho_beta) + (3 - 4 * kappa) * pow(rho_beta, 2) + pow(rho_beta, 3)) * exp(-rho_beta));
			}
			if (overlap == "1pz2pz") {
				Smatrix(k, l) = pow(3, 0.5) / (pow(1 - pow(tau, 2), 0.5) * tau * pow(rho, 3))
				                * (-pow(1 - kappa, 2) * (12 * (1 + kappa) * (1 + rho_alpha + 0.5 * pow(rho_alpha, 2)) + 2 * pow(rho_alpha, 3)) * exp(-rho_alpha)
				                   + (1 + kappa) * (12 * pow(1 - kappa, 2) * (1 + rho_beta + 0.5 * pow(rho_beta, 2)) + (3 - 4 * kappa) * pow(rho_beta, 3) + pow(rho_beta, 4)) * exp(-rho_beta));
			}
			if (overlap == "2pz2pz") {//wip
				Smatrix(k, l) = 1 / (pow(1 - pow(tau, 2), 0.5) * tau * pow(rho, 3))
				                * (-pow(1 - kappa, 2) * (48 * pow(1 + kappa, 2) * (1 + rho_alpha + 0.5 * pow(rho_alpha, 2)) + 2 * (5 + 6 * kappa) * pow(rho_alpha, 3) + 2 * pow(rho_alpha, 4)) * exp(-rho_alpha)
				                   + pow(1 + kappa, 2) * (48 * pow(1 - kappa, 2) * (1 + rho_beta + 0.5 * pow(rho_beta, 2)) + 2 * (5 - 6 * kappa) * pow(rho_beta, 3) + 2 * pow(rho_beta, 4)) * exp(-rho_beta));
			}
			if (overlap == "1pz2pz") {
				Smatrix(k, l) = (pow((1 - pow * (tau, 2)), 0.5) / (pow(6, 0.5) * (tau * rho))
				                 * (-(1 - kappa) * (1 - 2 * kappa) * exp(-rho_alpha) + ((1 - kappa) * (1 - 2 * kappa) + 2 * (1 - kappa) * (2 * (1 - kappa) + rho_beta)) * exp(-rho_beta)));
			}
			if (overlap == "2pz2pz") {
				Smatrix(k, l) = (pow((1 - pow * (tau, 2)), 0.5) / (pow(6, 0.5) * (tau * rho))
				                 * (-(1 - kappa) * (1 - 2 * kappa) * exp(-rho_alpha) + ((1 - kappa) * (1 - 2 * kappa) + 2 * (1 - kappa) * (2 * (1 - kappa) + rho_beta)) * exp(-rho_beta)));
			}

		}
	}
	return sum;

}

double distance (double x1, double y1, double z1, double x1, double y1, double z1) {
	return sqrt((pow((x1 - x2), 2) + pow((y1 - y2), 2) + pow((z1 - z2), 2)));
}
