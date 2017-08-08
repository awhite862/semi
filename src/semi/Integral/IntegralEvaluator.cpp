#   define ARMA_DONT_USE_WRAPPER
#include "IntegralEvaluator.h"
#include "semi/semi_utils.h"
#include <math.h>
using namespace arma;

namespace Semi {

//S_uv
void calculateOverlapMatrix(BasisSet<STOFunction> a, arma::mat &Smatrix) {
    Smatrix.set_size(a.myBasis.size(), a.myBasis.size());
    for (unsigned k = 0; k < a.myBasis.size(); k++) {
        for (unsigned l = 0; l < a.myBasis.size(); l++) {
            if (k == l) {
                Smatrix(k, l) =  1;
            }
            else if (a.myBasis[k].nlm.l + a.myBasis[l].nlm.l == 1) {
                if (a.myBasis[k].nlm.l == 0 && k > l) {
                    Smatrix(k, l) = -calculateOverlapSTO(a.myBasis[k], a.myBasis[l]);
                }
                else if (a.myBasis[k].nlm.l == 0 && k < l) {
                    Smatrix(k, l) = -calculateOverlapSTO(a.myBasis[k], a.myBasis[l]);
                }
                else {
                    Smatrix(k, l) = calculateOverlapSTO(a.myBasis[k], a.myBasis[l]);
                }
            }
            else {
                Smatrix(k, l) = calculateOverlapSTO(a.myBasis[k], a.myBasis[l]);
            }
        }
    }

    unsigned k = 0;
    while (k < a.myBasis.size()) {
        unsigned kl = a.myBasis[k].nlm.l;
        unsigned l = 0;
        while (l < a.myBasis.size()) {
            unsigned ll = a.myBasis[l].nlm.l;
            if (kl == 1 && ll == 1) {
                arma::mat rotationMatrix;
                findRotation(a.myBasis[k].x, a.myBasis[k].y, a.myBasis[k].z, a.myBasis[l].x, a.myBasis[l].y, a.myBasis[l].z, rotationMatrix);
                arma::mat sub(3, 3);
                sub = Smatrix.submat(k, l, k + 2, l + 2);
                sub = trans(rotationMatrix) * sub  * (rotationMatrix);
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        Smatrix(k + i, l + j) = sub(i, j);
                    }
                }
            }
            else if (kl == 1 && ll == 0) {
                arma::mat rotationMatrix;
                findRotation(a.myBasis[k].x, a.myBasis[k].y, a.myBasis[k].z, a.myBasis[l].x, a.myBasis[l].y, a.myBasis[l].z, rotationMatrix);
                arma::mat sub(3, 1);
                sub = Smatrix.submat(k, l, k + 2, l);
                sub = trans(rotationMatrix) * sub;
                for (int i = 0; i < 3; i++) {
                    Smatrix(k  + i, l) = sub(i, 0);
                }
            }
            else if (kl == 0 && ll == 1) {
                arma::mat rotationMatrix;
                findRotation(a.myBasis[k].x, a.myBasis[k].y, a.myBasis[k].z, a.myBasis[l].x, a.myBasis[l].y, a.myBasis[l].z, rotationMatrix);
                arma::mat sub(1, 3);
                sub = (Smatrix.submat(k, l, k, l + 2));
                sub = (sub * (rotationMatrix));
                for (int i = 0; i < 3; i++) {
                    Smatrix(k, l  + i) = sub(0, i);
                }
            }
            else {// do nothing for both s-functions
            }
            l += 2 * ll + 1;
        }
        k += 2 * kl + 1;
    }
}





double calculateOverlapCGTO(CGTOFunction a, CGTOFunction b) {
    double overlap = 0.0;
    for (int k = 0; k < 3; k++) {
        for (int i = 0; i < 3; i++) {
            GTOFunction tempA(a.nlm, a.a, a.b, a.c, a.alphaVec[k], a.r);
            GTOFunction tempB(b.nlm, b.a, b.b, b.c, b.alphaVec[i], b.r);
            overlap += a.nVec[k] * b.nVec[i] * calculateOverlapGTOUnnorm(tempA, tempB);
        }
    }
    return overlap;
}

double calculateOverlapGTO(GTOFunction a, GTOFunction b) {
    double a_i, b_j;
    if (a.a == 1) {
        a_i = 0;
    }
    else if (a.b == 1) {
        a_i = 1;
    }
    else if (a.c == 1) {
        a_i = 2;
    }
    if (b.a == 1) {
        b_j = 0;
    }
    else if (b.b == 1) {
        b_j = 1;
    }
    else if (b.c == 1) {
        b_j = 2;
    }
    if (a.nlm.l == 0 && b.nlm.l == 0) {
        return a.n * b.n * pow(M_PI / (a.alpha + b.alpha), 3.0 / 2.0) * exp((-a.alpha * b.alpha) / (a.alpha + b.alpha) * pow(distance(a.r(0), a.r(1), a.r(2), b.r(0), b.r(1), b.r(2)), 2));
    }
    else if (a.nlm.l == 1 && b.nlm.l == 0) {
        return -1.0 * b.alpha / (a.alpha + b.alpha) * (a.r(a_i) - b.r(a_i)) *  a.n * b.n * pow(M_PI / (a.alpha + b.alpha), 3.0 / 2.0) * exp((-a.alpha * b.alpha) / (a.alpha + b.alpha) * pow(distance(a.r(0), a.r(1), a.r(2), b.r(0), b.r(1), b.r(2)), 2));
    }
    else if (a.nlm.l == 0 && b.nlm.l == 1) {
        return -1.0 * a.alpha / (a.alpha + b.alpha) * (b.r(b_j) - a.r(b_j)) *  a.n * b.n * pow(M_PI / (a.alpha + b.alpha), 3.0 / 2.0) * exp((-a.alpha * b.alpha) / (a.alpha + b.alpha) * pow(distance(a.r(0), a.r(1), a.r(2), b.r(0), b.r(1), b.r(2)), 2));
    }
    else if (a.nlm.l == 1 && b.nlm.l == 1) {
        return (delta(a_i, b_j) / (2 * a.alpha + 2 * b.alpha) + (a.alpha * b.alpha) / pow(a.alpha + b.alpha, 2) * (a.r(a_i) - b.r(a_i)) * (b.r(b_j) - a.r(b_j))) * a.n * b.n * pow(M_PI / (a.alpha + b.alpha), 3.0 / 2.0) * exp((-a.alpha * b.alpha) / (a.alpha + b.alpha) * pow(distance(a.r(0), a.r(1), a.r(2), b.r(0), b.r(1), b.r(2)), 2));
    }
}

double calculateOverlapGTOUnnorm(GTOFunction a, GTOFunction b) {
    double a_i, b_j;
    if (a.a == 1) {
        a_i = 0;
    }
    else if (a.b == 1) {
        a_i = 1;
    }
    else if (a.c == 1) {
        a_i = 2;
    }
    if (b.a == 1) {
        b_j = 0;
    }
    else if (b.b == 1) {
        b_j = 1;
    }
    else if (b.c == 1) {
        b_j = 2;
    }
    if (a.nlm.l == 0 && b.nlm.l == 0) {
        return pow(M_PI / (a.alpha + b.alpha), 3.0 / 2.0) * exp((-a.alpha * b.alpha) / (a.alpha + b.alpha) * pow(distance(a.r(0), a.r(1), a.r(2), b.r(0), b.r(1), b.r(2)), 2));
    }
    else if (a.nlm.l == 1 && b.nlm.l == 0) {
        return -1.0 * b.alpha / (a.alpha + b.alpha) * (a.r(a_i) - b.r(a_i)) * pow(M_PI / (a.alpha + b.alpha), 3.0 / 2.0) * exp((-a.alpha * b.alpha) / (a.alpha + b.alpha) * pow(distance(a.r(0), a.r(1), a.r(2), b.r(0), b.r(1), b.r(2)), 2));
    }
    else if (a.nlm.l == 0 && b.nlm.l == 1) {
        return -1.0 * a.alpha / (a.alpha + b.alpha) * (b.r(b_j) - a.r(b_j)) * pow(M_PI / (a.alpha + b.alpha), 3.0 / 2.0) * exp((-a.alpha * b.alpha) / (a.alpha + b.alpha) * pow(distance(a.r(0), a.r(1), a.r(2), b.r(0), b.r(1), b.r(2)), 2));
    }
    else if (a.nlm.l == 1 && b.nlm.l == 1) {
        return (delta(a_i, b_j) / (2 * a.alpha + 2 * b.alpha) + (a.alpha * b.alpha) / pow(a.alpha + b.alpha, 2) * (a.r(a_i) - b.r(a_i)) * (b.r(b_j) - a.r(b_j))) * pow(M_PI / (a.alpha + b.alpha), 3.0 / 2.0) * exp((-a.alpha * b.alpha) / (a.alpha + b.alpha) * pow(distance(a.r(0), a.r(1), a.r(2), b.r(0), b.r(1), b.r(2)), 2));
    }
}

double calculateOverlapSTO(STOFunction a, STOFunction b) {
    double r = distance(a.x, a.y, a.z, b.x, b.y, b.z);
    double tau = (a.zeta - b.zeta) / (a.zeta + b.zeta);
    double rho = 0.5 * (a.zeta + b.zeta) * r;
    double kappa = 0.5 * (tau + 1.0 / tau);
    double rho_alpha = a.zeta * r, rho_beta = b.zeta * r;
    return calculateOverlapSTO(tau, rho, kappa, rho_alpha, rho_beta, a.nlm, b.nlm);
}

double calculateOverlapSTO(double tau, double rho, double kappa, double rho_alpha, double rho_beta, QNumber aQNumber, QNumber bQNumber) {
    if (fabs(rho) < tolerance) {
        return calculateOverlapSTOSamePosition(tau, rho, kappa, rho_alpha, rho_beta, aQNumber, bQNumber);
    } else if (fabs(tau) < tolerance) {
        return calculateOverlapSTOSameZeta(tau, rho, kappa, rho_alpha, rho_beta, aQNumber, bQNumber);
    }
    else {
        return calculateOverlapSTOFull(tau, rho, kappa, rho_alpha, rho_beta, aQNumber, bQNumber);
    }
}

double calculateOverlapSTOFull(double tau, double rho, double kappa, double rho_alpha, double rho_beta, QNumber aQNumber, QNumber bQNumber) {
    if (aQNumber.n == 1 && aQNumber.l == 0 && bQNumber.n == 1 && bQNumber.l == 0) { //1s1s s
        return (pow((1.0 - pow(tau, 2)), 0.5) / (tau * rho))
               * (-(1.0 - kappa) * (2.0 * (1.0 + kappa) + rho_alpha) * exp(-rho_alpha)
                  + (1.0 + kappa) * (2.0 * (1.0 - kappa) + rho_beta) * exp(-rho_beta));
    }
    if (aQNumber.n == 1 && aQNumber.l == 0 && bQNumber.n == 2 && bQNumber.l == 0) { //1s2s s
        return (pow((1.0 - pow(tau, 2)), 0.5) / (pow(3, 0.5) * tau * rho))
               * (-(1.0 - kappa) * (2.0 * (1.0 + kappa) * (2.0 - 3.0 * kappa) + (1.0 - 2.0 * kappa) * rho_alpha) * exp(-rho_alpha)
                  + (1.0 + kappa) * (2.0 * (1.0 - kappa) * (2.0 - 3.0 * kappa) + 4.0 * (1.0 - kappa) * rho_beta + pow(rho_beta, 2)) * exp(-rho_beta));
    }
    if (aQNumber.n == 2 && aQNumber.l == 0 && bQNumber.n == 2 && bQNumber.l == 0) { //2s2s s
        return (pow((1.0 - pow(tau, 2)), 0.5) / (3.0 * tau * rho))
               * (-(1.0 - kappa) * (2.0 * (1.0 + kappa) * (7 - 12.0 * pow(kappa, 2)) + 4.0 * (1.0 + kappa) * (2.0 - 3.0 * kappa) * rho_alpha + (1.0 - 2.0 * kappa) * pow(rho_alpha, 2)) *  exp(-rho_alpha)
                  + (1.0 + kappa) * (2.0 * (1.0 - kappa) * (7 - 12.0 * pow(kappa, 2)) + 4.0 * (1.0 - kappa) * (2.0 + 3.0 * kappa) * rho_beta + (1.0 + 2.0 * kappa) * pow(rho_beta, 2)) *  exp(-rho_beta));
    }
    if (aQNumber.n == 1 && aQNumber.l == 0 && bQNumber.n == 2 && bQNumber.l == 1 && aQNumber.m == bQNumber.m) { //1s2pz
        return pow((1.0 + tau) / (1.0 - tau), 0.5) * (1.0 / (tau * pow(rho, 2)))
               * (-pow(1.0 - kappa, 2) * (6.0 * (1.0 + kappa) * (1.0 + rho_alpha) + 2.0 * pow(rho_alpha, 2)) * exp(-rho_alpha)
                  + pow(1.0 + kappa, 1) * (6.0 * pow(1.0 - kappa, 2) * (1.0 + rho_beta) + 4.0 * (1.0 - kappa) * pow(rho_beta, 2) + pow(rho_beta, 3)) * exp(-rho_beta));
    }
    if (aQNumber.n == 2 && aQNumber.l == 0 && bQNumber.n == 2 && bQNumber.l == 1 && aQNumber.m == bQNumber.m) { //2s2pz s
        return pow((1.0 + tau) / (1.0 - tau), 0.5) * (1.0 / (pow(3, 0.5) * tau * pow(rho, 2)))
               * ((-pow(1.0 - kappa, 2) * (6.0 * (1.0 + kappa) * (3.0 + 4.0 * kappa) * (1.0 + rho_alpha) + 2.0 * (5.0 + 6.0 * kappa) * pow(rho_alpha, 2) + 2.0 * pow(rho_alpha, 3)) * exp(-rho_alpha))
                  + (1.0 + kappa) * (6.0 * pow(1.0 - kappa, 2) * (3.0 + 4.0 * kappa) * (1.0 + rho_beta) + 4.0 * (1.0 - kappa) * (2.0 + 3.0 * kappa) * pow(rho_beta, 2) + (1.0 + 2.0 * kappa) * pow(rho_beta, 3)) * exp(-rho_beta));
    }
    if (aQNumber.n == 2 && aQNumber.l == 1 && bQNumber.n == 2 && bQNumber.l == 1 && aQNumber.m == bQNumber.m && aQNumber.m == 0) { //2pz2pz s
        return -1.0 / (pow(1.0 - pow(tau, 2), 0.5) * tau * pow(rho, 3))
               * (-pow(1.0 - kappa, 2) * (48.0 * pow(1.0 + kappa, 2) * (1.0 + rho_alpha + 0.5 * pow(rho_alpha, 2)) + 2.0 * (5.0 + 6.0 * kappa) * pow(rho_alpha, 3) + 2.0 * pow(rho_alpha, 4)) * exp(-rho_alpha)
                  + pow(1.0 + kappa, 2) * (48.0 * pow(1.0 - kappa, 2) * (1.0 + rho_beta + 0.5 * pow(rho_beta, 2)) + 2.0 * (5.0 - 6.0 * kappa) * pow(rho_beta, 3) + 2.0 * pow(rho_beta, 4)) * exp(-rho_beta));
    }
    if (aQNumber.n == 2 && aQNumber.l == 1 && bQNumber.n == 2 && bQNumber.l == 1 && aQNumber.m == bQNumber.m) { //2pz2pz p
        return 1.0 / (pow(1.0 - pow(tau, 2), 0.5) * tau * pow(rho, 3))
               * (-pow(1.0 - kappa, 2) * (24.0 * pow(1.0 + kappa, 2) * (1.0 + rho_alpha) + 12.0 * (1.0 + kappa) * pow(rho_alpha, 2) + 2 * pow(rho_alpha, 3)) * exp(-rho_alpha)
                  + pow(1.0 + kappa, 2) * (24.0 * pow(1.0 - kappa, 2) * (1.0 + rho_beta) + 12.0 * (1.0 - kappa) * pow(rho_beta, 2) + 2 * pow(rho_beta, 3)) * exp(-rho_beta));
    }
    if (isReversed(aQNumber, bQNumber)) {
        return calculateOverlapSTOFull(-tau, rho, -kappa, rho_beta, rho_alpha, bQNumber, aQNumber);
    }
    return 0;
}

double calculateOverlapSTOSameZeta(double tau, double rho, double kappa, double rho_alpha, double rho_beta, QNumber aQNumber, QNumber bQNumber) { //Tau = 0
    if (aQNumber.n == 1 && aQNumber.l == 0 && bQNumber.n == 1 && bQNumber.l == 0) { //1s1s
        return (1.0 + rho + 1.0 / 3.0 * pow(rho, 2)) * exp(-rho);
    }
    if (aQNumber.n == 1 && aQNumber.l == 0 && bQNumber.n == 2 && bQNumber.l == 0) { //1s2s
        return pow(3, 0.5) / 2.0 * (1.0 + rho + 4 / 9 * pow(rho, 2) + 1.0 / 9 * pow(rho, 3)) * exp(-rho);
    }
    if (aQNumber.n == 2 && aQNumber.l == 0 && bQNumber.n == 2 && bQNumber.l == 0) { //2s2s
        return (1.0 + rho + 4.0 / 9.0 * pow(rho, 2) + 1.0 / 9.0 * pow(rho, 3) + 1.0 / 45 * pow(rho, 4)) * exp(-rho);
    }
    if (aQNumber.n == 1 && aQNumber.l == 0 && bQNumber.n == 2 && bQNumber.l == 1 && aQNumber.m == bQNumber.m) { //1s2pz
        return  1.0 / 2.0 * rho * (1.0 + rho + 1.0 / 3.0 * pow(rho, 2)) * exp(-rho);
    }
    if (aQNumber.n == 2 && aQNumber.l == 0 && bQNumber.n == 2 && bQNumber.l == 1 && aQNumber.m == bQNumber.m) { //2s2pz
        return  1.0 / (2.0 * pow(3, 0.5)) * rho * (1.0 + rho + 7.0 / 15.0 * pow(rho, 2) + 2.0 / 15.0 * pow(rho, 3)) * exp(-rho);
    }
    if (aQNumber.n == 2 && aQNumber.l == 1 && bQNumber.n == 2 && bQNumber.l == 1 && aQNumber.m == bQNumber.m && aQNumber.m == 0) { //2pz2pz  
        return  -(-1.0 - rho - 1.0 / 5.0 * pow(rho, 2) + 2.0 / 15.0 * pow(rho, 3) + 1.0 / 15.0 * pow(rho, 4)) * exp(-rho);
    }
    if (aQNumber.n == 2 && aQNumber.l == 1 && bQNumber.n == 2 && bQNumber.l == 1 && aQNumber.m == bQNumber.m) { //2pz2pz
        return  (1.0 + rho + 2.0 / 5.0 * pow(rho, 2) + 1.0 / 15.0 * pow(rho, 3)) * exp(-rho);
    }
    if (isReversed(aQNumber, bQNumber)) {
        return calculateOverlapSTOSameZeta(-tau, rho, -kappa, rho_beta, rho_alpha, bQNumber, aQNumber);
    }
    return 0;
}

double calculateOverlapSTOSamePosition(double tau, double rho, double kappa, double rho_alpha, double rho_beta, QNumber aQNumber, QNumber bQNumber) { //Rho = 0
    if (aQNumber.n == 1 && aQNumber.l == 0 && bQNumber.n == 1 && bQNumber.l == 0) { //1s1s
        return pow(1.0 + tau, 1.5) * pow(1.0 - tau, 1.5);
    }
    if (aQNumber.n == 1 && aQNumber.l == 0 && bQNumber.n == 2 && bQNumber.l == 0) { //1s2s
        return 1.0 / 2.0 * pow(3, 0.5) * pow(1.0 + tau, 1.5) * pow(1.0 - tau, 2.5);
    }
    if (aQNumber.n == 2 && aQNumber.l == 0 && bQNumber.n == 2 && bQNumber.l == 0) { //2s2s
        return pow(1.0 + tau, 2.5) * pow(1.0 - tau, 2.5);
    }
    if (aQNumber.n == 1 && aQNumber.l == 0 && bQNumber.n == 2 && bQNumber.l == 1 && aQNumber.m == bQNumber.m) { //1s2pz
        return 0;
    }
    if (aQNumber.n == 2 && aQNumber.l == 0 && bQNumber.n == 2 && bQNumber.l == 1 && aQNumber.m == bQNumber.m) { //2s2pz
        return 0;
    }
    if (aQNumber.n == 2 && aQNumber.l == 1 && bQNumber.n == 2 && bQNumber.l == 1 && aQNumber.m == bQNumber.m && aQNumber.m == 0) { //2pz2pz s
        return pow(1.0 + tau, 2.5) * pow(1.0 - tau, 2.5);
    }
    if (aQNumber.n == 2 && aQNumber.l == 1 && bQNumber.n == 2 && bQNumber.l == 1 && aQNumber.m == bQNumber.m) { //2pz2pz p
        return pow(1.0 + tau, 2.5) * pow(1.0 - tau, 2.5);
    }
    if (isReversed(aQNumber, bQNumber)) {
        return calculateOverlapSTOSamePosition(-tau, rho, -kappa, rho_beta, rho_alpha, bQNumber, aQNumber);
    }
    return 0;
}

double calculateBasicCoulombIntegral(double zeta, double tau, double rho, double kappa, double rho_alpha, double rho_beta, QNumber aQNumber, QNumber bQNumber) {
    if (std::abs(rho) < tolerance) {
        return calculateBasicCoulombIntegralSamePosition(zeta, tau, rho, kappa, rho_alpha, rho_beta, aQNumber, bQNumber);
    } else if (std::abs(tau) < tolerance) {
        return calculateBasicCoulombIntegralSameZeta(zeta, tau, rho, kappa, rho_alpha, rho_beta, aQNumber, bQNumber);
    }
    else {
        return calculateBasicCoulombIntegralFull(zeta, tau, rho, kappa, rho_alpha, rho_beta, aQNumber, bQNumber);
    }
}

double calculateBasicCoulombIntegralFull(double zeta, double tau, double rho, double kappa, double rho_alpha, double rho_beta, QNumber aQNumber, QNumber bQNumber) {
    if (aQNumber.n == 1 && aQNumber.l == 0 && bQNumber.n == 1 && bQNumber.l == 0) { //1s1s
        return (zeta / rho) * (1.0 - pow(1 - kappa, 2.0) * (1.0 / 4.0 * (2.0 + kappa) + 1.0 / 4.0 * rho_alpha) * exp(-2.0 * rho_alpha)
                                 - (1.0 + pow(kappa, 2.0)) * (1.0 / 4.0 * (2.0 - kappa) + 1.0 / 4.0 * rho_beta) * exp(-2.0 * rho_beta));
    }
    if (aQNumber.n == 1 && aQNumber.l == 0 && bQNumber.n == 2 && bQNumber.l == 0) { //1s2s
        return (zeta / rho) * (1.0 - (pow(1.0 - kappa, 3.0) * (1.0 / 16.0 * (1.0 - 5.0 * kappa - 4.0 * pow(kappa, 2.0)) - 1.0 / 8.0 * kappa * rho_alpha) * exp(-2.0 * rho_alpha))
                                 - (pow(1.0 + kappa, 2.0) * (1.0 / 16.0 * (15.0 - 22.0 * kappa + 15.0 * pow(kappa, 2.0) - 4.0 * pow(kappa, 3.0))
                                         + 3.0 / 8.0 * (3.0 - 3.0 * kappa + pow(kappa, 2.0)) * rho_beta + 1.0 / 4.0 * (2.0 - kappa) * pow(rho_beta, 2.0) + 1.0 / 12.0 * pow(rho_beta, 3.0))) * exp(-2.0 * rho_beta));
    }
    if (aQNumber.n == 2 && aQNumber.l == 0 && bQNumber.n == 1 && bQNumber.l == 0) { //1s2s
        kappa = -kappa;
        return (zeta / rho) * (1.0 - (pow(1.0 - kappa, 3.0) * (1.0 / 16.0 * (1.0 - 5.0 * kappa - 4.0 * pow(kappa, 2.0)) - 1.0 / 8.0 * kappa * rho_beta) * exp(-2.0 * rho_beta))
                                 - (pow(1.0 + kappa, 2.0) * (1.0 / 16.0 * (15.0 - 22.0 * kappa + 15.0 * pow(kappa, 2.0) - 4.0 * pow(kappa, 3.0))
                                         + 3.0 / 8.0 * (3.0 - 3.0 * kappa + pow(kappa, 2.0)) * rho_alpha + 1.0 / 4.0 * (2.0 - kappa) * pow(rho_alpha, 2.0) + 1.0 / 12.0 * pow(rho_alpha, 3.0))) * exp(-2.0 * rho_alpha));
    }
    if (aQNumber.n == 2 && aQNumber.l == 0 && bQNumber.n == 2 && bQNumber.l == 0) { //2s2s
        return (zeta / rho) * (1.0 - pow(1 - kappa, 3.0) * (1.0 / 16.0 * (8.0 - kappa - 27.0 * pow(kappa, 2.0) - 30.0 * pow(kappa, 3.0) - 10.0 * pow(kappa, 4.0))
                                 + 1.0 / 32.0 * (11.0 - 19.0 * kappa - 44.0 * pow(kappa, 2) - 20.0 * pow(kappa, 3)) * rho_alpha
                                 + 1.0 / 16.0 * (1.0 - 5.0 * kappa - 4.0 * pow(kappa, 2.0)) * pow(rho_alpha, 2) - 1.0 / 24.0 * kappa * pow(rho_alpha, 3)) * exp(-2.0 * rho_alpha)
                                 - pow(1.0 + kappa, 3.0) * (1.0 / 16.0 * (8.0 + kappa - 27.0 * pow(kappa, 2.0) + 30.0 * pow(kappa, 3.0) - 10.0 * pow(kappa, 4))
                                         + 1.0 / 32.0 * (11.0 + 19.0 * kappa - 44.0 * pow(kappa, 2.0) + 20.0 * pow(kappa, 3)) * rho_beta
                                         + 1.0 / 16.0 * (1.0 + 5.0 * kappa - 4.0 * pow(kappa, 2.0)) * pow(rho_beta, 2) + 1.0 / 24.0 * kappa * pow(rho_beta, 3)) * exp(-2 * rho_beta));
    }
    return 0;
}

//[1Sa,1Sb] tau = 0
double calculateBasicCoulombIntegralSameZeta(double zeta, double tau, double rho, double kappa, double rho_alpha, double rho_beta, QNumber aQNumber, QNumber bQNumber) {
    if (aQNumber.n == 1 && aQNumber.l == 0 && bQNumber.n == 1 && bQNumber.l == 0) { //1s1s
        return (zeta / rho) * (1.0 - (1.0 + 11.0 / 8.0 * rho + 3.0 / 4.0 * pow(rho, 2.0) + 1.0 / 6.0 * pow(rho, 3.0)) * exp(-2.0 * rho));
    }
    if (aQNumber.n == 1 && aQNumber.l == 0 && bQNumber.n == 2 && bQNumber.l == 0) { //1s2s
        return (zeta / rho) * (1.0 - (1.0 + 25.0 / 16.0 * rho + 9.0 / 8.0 * pow(rho, 2.0) + 23.0 / 48.0 * pow(rho, 3.0) + 1.0 / 8.0 * pow(rho, 4.0) + 1.0 / 60.0 * pow(rho, 5.0)) * exp(-2.0 * rho));
    }
    if (aQNumber.n == 2 && aQNumber.l == 0 && bQNumber.n == 1 && bQNumber.l == 0) { //1s2s
        return (zeta / rho) * (1.0 - (1.0 + 25.0 / 16.0 * rho + 9.0 / 8.0 * pow(rho, 2.0) + 23.0 / 48.0 * pow(rho, 3.0) + 1.0 / 8.0 * pow(rho, 4.0) + 1.0 / 60.0 * pow(rho, 5.0)) * exp(-2.0 * rho));
    }
    if (aQNumber.n == 2 && aQNumber.l == 0 && bQNumber.n == 2 && bQNumber.l == 0) { //2s2s
        return (zeta / rho) * (1.0 - (1.0 + 419.0 / 256.0 * rho + 163.0 / 128.0 * pow(rho, 2.0) + 119.0 / 192.0 * pow(rho, 3.0) + 5.0 / 24.0 * pow(rho, 4) + 1.0 / 20.0 * pow(rho, 5) + 1.0 / 120.0 * pow(rho, 6) + 1.0 / 1260.0 * pow(rho, 7)) * exp(-2.0 * rho));
    }
    return 0;
}

//[1Sa,1Sb] rho = 0
double calculateBasicCoulombIntegralSamePosition(double zeta, double tau, double rho, double kappa, double rho_alpha, double rho_beta, QNumber aQNumber, QNumber bQNumber) {
    if (aQNumber.n == 1 && aQNumber.l == 0 && bQNumber.n == 1 && bQNumber.l == 0) { //1s1s
        return 1.0 / 8.0 * (1.0 - pow(tau, 2.0)) * (5.0 - pow(tau, 2.0)) * zeta;
    }
    if (aQNumber.n == 1 && aQNumber.l == 0 && bQNumber.n == 2 && bQNumber.l == 0) { //1s2s
        return 1.0 / 32.0 * (1.0 - pow(tau, 2.0)) * (14.0 - 7 * tau - 3.0 * pow(tau, 2.0) + pow(tau, 3)) * zeta;
    }
    if (aQNumber.n == 2 && aQNumber.l == 0 && bQNumber.n == 1 && bQNumber.l == 0) { //1s2s
        return 1.0 / 32.0 * (1.0 - pow(-tau, 2.0)) * (14.0 - 7 * -tau - 3.0 * pow(-tau, 2.0) + pow(-tau, 3)) * zeta;
    }
    if (aQNumber.n == 2 && aQNumber.l == 0 && bQNumber.n == 2 && bQNumber.l == 0) { //2s2s
        return 1.0 / 256.0 * (1.0 - pow(tau, 2.0)) * (93.0 - 47.0 * pow(tau, 2.0) + 23.0 * pow(tau, 4) - 5.0 * pow(tau, 6)) * zeta;
    }
    return 0;
}

void findRotation(double x1, double y1, double z1, double x2, double y2, double z2, arma::mat &rotationMatrix) {
    double x = x2 - x1;
    double y = y2 - y1;
    double z = z2 - z1;

    arma::mat R1(3, 3); R1.eye();
    arma::mat R2(3, 3); R2.eye();

    // rotate about z if y != 0
    if (std::abs(y) > tolerance) {
        double theta = std::atan2(y, x);
        double c = std::cos(theta);
        double s = std::sin(theta);
        R1(0, 0) = c;
        R1(1, 1) = c;
        R1(0, 1) = s;
        R1(1, 0) = -s;
    }

    arma::vec v(3); v(0) = x; v(1) = y; v(2) = z;
    arma::vec u = R1 * v;

    // finally rotate about y
    double theta = std::atan2(u(0), u(2));
    double c1 = std::cos(theta);
    double s1 = std::sin(theta);
    R2(0, 0) = c1;
    R2(2, 2) = c1;
    R2(0, 2) = -s1;
    R2(2, 0) = s1;

    rotationMatrix = R2 * R1;
}

bool isReversed(QNumber aQNumber, QNumber bQNumber) {
    if ((bQNumber.n == 0 && bQNumber.l == 0 && aQNumber.n == 1 && aQNumber.l == 0) ||
            (bQNumber.n == 1 && bQNumber.l == 0 && aQNumber.n == 1 && aQNumber.l == 0) ||
            (bQNumber.n == 0 && bQNumber.l == 0 && aQNumber.n == 2 && aQNumber.l == 0) ||
            (bQNumber.n == 1 && bQNumber.l == 0 && aQNumber.n == 2 && aQNumber.l == 0) ||
            (bQNumber.n == 2 && bQNumber.l == 0 && aQNumber.n == 2 && aQNumber.l == 0) ||
            (bQNumber.n == 0 && bQNumber.l == 0 && aQNumber.n == 2 && aQNumber.l == 1 && aQNumber.m == bQNumber.m) ||
            (bQNumber.n == 1 && bQNumber.l == 0 && aQNumber.n == 2 && aQNumber.l == 1 && aQNumber.m == bQNumber.m) ||
            (bQNumber.n == 2 && bQNumber.l == 0 && aQNumber.n == 2 && aQNumber.l == 1 && aQNumber.m == bQNumber.m) ||
            (bQNumber.n == 1 && bQNumber.l == 0 && aQNumber.n == 1 && aQNumber.l == 1 && aQNumber.m == bQNumber.m) ||
            (bQNumber.n == 1 && bQNumber.l == 1 && aQNumber.n == 2 && aQNumber.l == 0 && aQNumber.m == bQNumber.m) ||
            (bQNumber.n == 1 && bQNumber.l == 1 && aQNumber.n == 2 && aQNumber.l == 1 && aQNumber.m == bQNumber.m) ||
            (bQNumber.n == 2 && bQNumber.l == 1 && aQNumber.n == 2 && aQNumber.l == 1 && aQNumber.m == bQNumber.m) ||
            (bQNumber.n == 1 && bQNumber.l == 1 && aQNumber.n == 2 && aQNumber.l == 1 && aQNumber.m == bQNumber.m) ||
            (bQNumber.n == 2 && bQNumber.l == 1 && aQNumber.n == 2 && aQNumber.l == 1 && aQNumber.m == bQNumber.m)) {
        return true;
    }
    return false;
}


} //namespace Semi
