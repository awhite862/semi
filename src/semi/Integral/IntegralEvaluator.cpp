#define ARMA_DONT_USE_WRAPPER
#include "IntegralEvaluator.h"
#include "semi/semi_utils.h"
#include <math.h>
using namespace arma;

namespace Semi {

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

double calculateOverlapSTOFull(double tau, double rho, double kappa, double rho_alpha, double rho_beta, int *a, int *b) {
    if (a[0] == 0 && a[1] == 0 && b[0] == 1 && b[1] == 0) { //0s1s
        return (pow((1.0 - pow(tau, 2)), 0.5) / (pow(2, 0.5) * tau * rho))
               * (-(1.0 - kappa) * exp(-rho_alpha) + ((1.0 - kappa) + rho_beta) * exp(-rho_beta));
    }
    if (a[0] == 1 && a[1] == 0 && b[0] == 1 && b[1] == 0) { //1s1s
        return (pow((1.0 - pow(tau, 2)), 0.5) / (tau * rho))
               * (-(1.0 - kappa) * (2.0 * (1.0 + kappa) + rho_alpha) * exp(-rho_alpha)
                  + (1.0 + kappa) * (2.0 * (1.0 - kappa) + rho_beta) * exp(-rho_beta));
    }
    if (a[0] == 0 && a[1] == 0 && b[0] == 2 && b[1] == 0) { //0s2s
        return (pow((1.0 - pow(tau, 2)), 0.5) / (pow(6, 0.5) * tau * rho))
               * (-(1.0 - kappa) * (1.0 - 2.0 * kappa) * exp(-rho_alpha)
                  + ((1.0 - kappa) * (1.0 - 2.0 * kappa) + 2.0 * (1.0 - kappa) * rho_beta + pow(rho_beta, 2)) * exp(-rho_beta));
    }
    if (a[0] == 1 && a[1] == 0 && b[0] == 2 && b[1] == 0) { //1s2s
        return (pow((1.0 - pow(tau, 2)), 0.5) / (pow(3, 0.5) * tau * rho))
               * (-(1.0 - kappa) * (2.0 * (1.0 + kappa) * (2.0 - 3.0 * kappa) + (1.0 - 2.0 * kappa) * exp(-rho_alpha))
                  + (1.0 + kappa) * (2.0 * (1.0 - kappa) * (2.0 - 3.0 * kappa) + 4.0 * (1.0 - kappa) * rho_beta + pow(rho_beta, 2)) * exp(-rho_beta));
    }
    if (a[0] == 2 && a[1] == 0 && b[0] == 2 && b[1] == 0) { //2s2s
        return (pow((1.0 - pow(tau, 2)), 0.5) / (3.0 * tau * rho))
               * (-(1.0 - kappa) * (2.0 * (1.0 + kappa) * (7 - 12.0 * pow(kappa, 2)) + 4.0 * (1.0 + kappa) * (2.0 - 3.0 * kappa) * rho_alpha + (1.0 - 2.0 * kappa) * pow(rho_alpha, 2) *  exp(-rho_alpha))
                  + (1.0 + kappa) * (2.0 * (1.0 - kappa) * (7 - 12.0 * pow(kappa, 2)) + 4.0 * (1.0 - kappa) * (2.0 + 3.0 * kappa) * rho_beta + (1.0 + 2.0 * kappa) * pow(rho_beta, 2) *  exp(-rho_beta)));
    }
    if (a[0] == 0 && a[1] == 0 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //0s2pz
        return pow((1.0 + tau) / (1.0 - tau), 0.5) * (1.0 / (pow(2, 0.5) * tau * pow(rho, 2)))
               * (-2.0 * pow(1.0 - kappa, 2) * (1.0 + rho_alpha) * exp(-rho_alpha)
                  + (2.0 * pow(1.0 - kappa, 2) * (1.0 + rho_beta) + 2.0 * (1.0 - kappa) * pow(rho_beta, 2) + pow(rho_beta, 3)) * exp(-rho_beta));
    }
    if (a[0] == 1 && a[1] == 0 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //1s2pz
        return pow((1.0 + tau) / (1.0 - tau), 0.5) * (1.0 / (tau * pow(rho, 2)))
               * (-pow(1.0 - kappa, 2) * (6.0 * (1.0 + kappa) * (1.0 + rho_alpha) + 2.0 * pow(rho_alpha, 2)) * exp(-rho_alpha)
                  + pow(1.0 + kappa, 2) * (6.0 * (1.0 - kappa) * (1.0 + rho_beta) + 4.0 * (1.0 - kappa) * pow(rho_beta, 2) + pow(rho_beta, 3)) * exp(-rho_beta));
    }
    if (a[0] == 2 && a[1] == 0 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //2s2pz
        return pow((1.0 + tau) / (1.0 - tau), 0.5) * (1.0 / (pow(3, 0.5) * tau * pow(rho, 2)))
               * (-pow(1.0 - kappa, 2) * (6.0 * (1.0 + kappa) * (3.0 + 4.0 * kappa) * (1.0 + rho_alpha) + 2.0 * (5.0 + 6.0 * kappa) * pow(rho_alpha, 2) + 2.0 * pow(rho_alpha, 3)) * exp(-rho_alpha)
                  + (1.0 + kappa) * (6.0 * pow(1.0 - kappa, 2) * (3.0 + 4.0 * kappa) * (1.0 + rho_beta) + 4.0 * (1.0 - kappa) * (2.0 + 3.0 * kappa) * pow(rho_beta, 2) + (1.0 + 2.0 * kappa) * pow(rho_alpha, 3)) * exp(-rho_beta));
    }
    if (a[0] == 1 && a[1] == 0 && b[0] == 1 && b[1] == 1 && a[2] == b[2]) { //1pz1s
        return pow((1.0 - tau) / (1.0 + tau), 0.5) * (pow(3, 0.5) / (tau * pow(rho, 2)))
               * (-(1.0 - kappa) * (2.0 * (1.0 + kappa) * (1.0 + rho_alpha) + pow(rho_alpha, 2)) * exp(-rho_alpha)
                  + (1.0 + kappa) * (2.0 * (1.0 - kappa) * (1.0 + rho_beta) + pow(rho_beta, 2)) * exp(-rho_beta));
    }
    if (a[0] == 1 && a[1] == 1 && b[0] == 2 && b[1] == 0 && a[2] == b[2]) { //1pz2s
        return pow((1.0 - tau) / (1.0 + tau), 0.5) * (1.0 / (tau * pow(rho, 2)))
               * (-(1.0 - kappa) * (2.0 * (1.0 + kappa) * (2.0 - 3.0 * kappa) * (1.0 + rho_alpha) + (1.0 - 2.0 * kappa) * pow(rho_alpha, 2)) * exp(-rho_alpha)
                  + (1.0 + kappa) * (2.0 * (1.0 - kappa) * (2.0 - 3.0 * kappa) * (1.0 + rho_beta) + (3.0 - 4.0 * kappa) * pow(rho_beta, 2) + pow(rho_beta, 3)) * exp(-rho_beta));
    }
    if (a[0] == 1 && a[1] == 1 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //1pz2pz
        return pow(3, 0.5) / (pow(1.0 - pow(tau, 2), 0.5) * tau * pow(rho, 3))
               * (-pow(1.0 - kappa, 2) * (12.0 * (1.0 + kappa) * (1.0 + rho_alpha + 0.5 * pow(rho_alpha, 2)) + 2.0 * pow(rho_alpha, 3)) * exp(-rho_alpha)
                  + (1.0 + kappa) * (12.0 * pow(1.0 - kappa, 2) * (1.0 + rho_beta + 0.5 * pow(rho_beta, 2)) + (3.0 - 4.0 * kappa) * pow(rho_beta, 3) + pow(rho_beta, 4)) * exp(-rho_beta));
    }
    if (a[0] == 2 && a[1] == 1 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //2pz2pz
        return 1.0 / (pow(1.0 - pow(tau, 2), 0.5) * tau * pow(rho, 3))
               * (-pow(1.0 - kappa, 2) * (48.0 * pow(1.0 + kappa, 2) * (1.0 + rho_alpha + 0.5 * pow(rho_alpha, 2)) + 2.0 * (5.0 + 6.0 * kappa) * pow(rho_alpha, 3) + 2.0 * pow(rho_alpha, 4)) * exp(-rho_alpha)
                  + pow(1.0 + kappa, 2) * (48.0 * pow(1.0 - kappa, 2) * (1.0 + rho_beta + 0.5 * pow(rho_beta, 2)) + 2.0 * (5.0 - 6.0 * kappa) * pow(rho_beta, 3) + 2.0 * pow(rho_beta, 4)) * exp(-rho_beta));
    }
    if (a[0] == 1 && a[1] == 1 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //1pz2pz
        return pow(3, 0.5) / (pow(1.0 - pow(tau, 2), 0.5) * tau * pow(rho, 3))
               * (-pow(1.0 - kappa, 2) * (6.0 * (1.0 + kappa) * (1.0 + rho_alpha) + 2.0 * pow(rho_alpha, 2)) * exp(-rho_alpha)
                  + (1.0 + kappa) * (6.0 * pow(1.0 - kappa, 2) * (1.0 + rho_beta) + 4.0 * (1.0 - kappa) * pow(rho_beta, 2) + pow(rho_beta, 3)) * exp(-rho_beta));
    }
    if (a[0] == 2 && a[1] == 1 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //2pz2pz
        return 1.0 / (pow(1.0 - pow(tau, 2), 0.5) * tau * pow(rho, 3))
               * (-pow(1.0 - kappa, 2) * (24.0 * pow(1.0 + kappa, 2) * (1.0 + rho_alpha) + 12.0 * (1.0 + kappa) * pow(rho_alpha, 2) + pow(rho_alpha, 3)) * exp(-rho_alpha)
                  + pow(1.0 + kappa, 2) * (24.0 * pow(1.0 - kappa, 2) * (1.0 + rho_beta) + 12.0 * (1.0 - kappa) * pow(rho_beta, 2) + pow(rho_beta, 3)) * exp(-rho_beta));
    }
    if (isReversed(a, b)) {
        return calculateOverlapSTOFull(-tau, rho, -kappa, rho_beta, rho_alpha, b, a);
    }
    return 0;
    //throw std::runtime_error("Failed to catch case");
}

double calculateOverlapSTOSameZeta(double tau, double rho, double kappa, double rho_alpha, double rho_beta, int *a, int *b) { //Tau = 0
    if (a[0] == 0 && a[1] == 0 && b[0] == 1 && b[1] == 0) { //0s1s
        return (1.0 / pow(2, 0.5)) * (1.0 + rho) * exp(-rho);
    }
    if (a[0] == 1 && a[1] == 0 && b[0] == 1 && b[1] == 0) { //1s1s
        return (1.0 + rho + 1.0 / 3.0 * pow(rho, 2)) * exp(-rho);
    }
    if (a[0] == 0 && a[1] == 0 && b[0] == 2 && b[1] == 0) { //0s2s
        return pow(6, -0.5) * (1.0 + rho + 2.0 / 3.0 * pow(rho, 2)) * exp(-rho);
    }
    if (a[0] == 1 && a[1] == 0 && b[0] == 2 && b[1] == 0) { //1s2s
        return pow(3, 0.5) / 2.0 * (1.0 + rho + 4 / 9 * pow(rho, 2) + 1.0 / 9 * pow(rho, 3)) * exp(-rho);
    }
    if (a[0] == 2 && a[1] == 0 && b[0] == 2 && b[1] == 0) { //2s2s
        return (1.0 + rho + 4 / 9 * pow(rho, 2) + 1.0 / 9 * pow(rho, 3) + 1.0 / 45 * pow(rho, 4)) * exp(-rho);
    }
    if (a[0] == 0 && a[1] == 0 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //0s2pz
        return  (1.0 / 3.0 * pow(2, 0.5) * rho) * (1.0 + rho) * exp(-rho);
    }
    if (a[0] == 1 && a[1] == 0 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //1s2pz
        return  1.0 / 2.0 * rho * (1.0 + rho + 1.0 / 3.0 * pow(rho, 2)) * exp(-rho);
    }
    if (a[0] == 2 && a[1] == 0 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //2s2pz
        return  1.0 / (2.0 * pow(3, 0.5)) * rho * (1.0 + rho + 7 / 15 * pow(rho, 2) + 2.0 / 15 * pow(rho, 3)) * exp(-rho);
    }
    if (a[0] == 1 && a[1] == 0 && b[0] == 1 && b[1] == 1 && a[2] == b[2]) { //1pz1s
        return  1.0 / pow(3, 0.5) * rho * (1.0 + rho) * exp(-rho);
    }
    if (a[0] == 1 && a[1] == 1 && b[0] == 2 && b[1] == 0 && a[2] == b[2]) { //1pz2s
        return  1.0 / 6.0 * (1.0 + rho + pow(rho, 2)) * exp(-rho);
    }
    if (a[0] == 1 && a[1] == 1 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //1pz2pz
        return  1.0 / 2.0 * pow(3, 0.5) * (-1.0 - rho + 1.0 / 3.0 * pow(rho, 3)) * exp(-rho);
    }
    if (a[0] == 2 && a[1] == 1 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //2pz2pz
        return  (-1.0 - rho - 1.0 / 5 * pow(rho, 2) + 2.0 / 15 * pow(rho, 3) + 1.0 / 15 * pow(rho, 4)) * exp(-rho);
    }
    if (a[0] == 1 && a[1] == 1 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //1pz2pz
        return 1.0 / 2.0 *  pow(3, 0.5) * (1.0 + rho + 1.0 / 3.0 * pow(rho, 2)) * exp(-rho);
    }
    if (a[0] == 2 && a[1] == 1 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //2pz2pz
        return  (1.0 + rho + 2.0 / 5 * pow(rho, 2) + 1.0 / 15 * pow(rho, 3)) * exp(-rho);
    }
    if (isReversed(a, b)) {
        return calculateOverlapSTOSameZeta(-tau, rho, -kappa, rho_alpha, rho_beta, b, a);
    }
    return 0;
    //throw std::runtime_error("Failed to catch case");
}

double calculateOverlapSTOSamePosition(double tau, double rho, double kappa, double rho_alpha, double rho_beta, int *a, int *b) { //Rho = 0
    if (a[0] == 0 && a[1] == 0 && b[0] == 1 && b[1] == 0) { //0s1s
        return 1.0 / pow(2, 0.5) * pow(1.0 + tau, 0.5) * pow(1.0 - tau, 1.5);
    }
    if (a[0] == 1 && a[1] == 0 && b[0] == 1 && b[1] == 0) { //1s1s
        return pow(1.0 + tau, 1.5) * pow(1.0 - tau, 1.5);
    }
    if (a[0] == 0 && a[1] == 0 && b[0] == 2 && b[1] == 0) { //0s2s
        return pow(6, -0.5) * pow(1.0 + tau, 0.5) * pow(1.0 - tau, 2.5);
    }
    if (a[0] == 1 && a[1] == 0 && b[0] == 2 && b[1] == 0) { //1s2s
        return 1.0 / 2.0 * pow(3, 0.5) * pow(1.0 + tau, 1.5) * pow(1.0 - tau, 2.5);
    }
    if (a[0] == 2 && a[1] == 0 && b[0] == 2 && b[1] == 0) { //2s2s
        return pow(1.0 + tau, 2.5) * pow(1.0 - tau, 2.5);
    }
    if (a[0] == 0 && a[1] == 0 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //0s2pz
        return 0;
    }
    if (a[0] == 1 && a[1] == 0 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //1s2pz
        return 0;
    }
    if (a[0] == 2 && a[1] == 0 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //2s2pz
        return 0;
    }
    if (a[0] == 1 && a[1] == 0 && b[0] == 1 && b[1] == 1 && a[2] == b[2]) { //1pz1s
        return 0;
    }
    if (a[0] == 1 && a[1] == 1 && b[0] == 2 && b[1] == 0 && a[2] == b[2]) { //1pz2s
        return 0;
    }
    if (a[0] == 1 && a[1] == 1 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //1pz2pz
        return 1.0 / 2.0 * pow(3, 0.5) * pow(1.0 + tau, 1.5) * pow(1.0 - tau, 2.5);
    }
    if (a[0] == 2 && a[1] == 1 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //2pz2pz
        return 1.0 / 2.0 * pow(3, 0.5) * pow(1.0 + tau, 1.5) * pow(1.0 - tau, 2.5);
    }
    if (a[0] == 1 && a[1] == 1 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //1pz2pz
        return pow(1.0 + tau, 2.5) * pow(1.0 - tau, 2.5);
    }
    if (a[0] == 2 && a[1] == 1 && b[0] == 2 && b[1] == 1 && a[2] == b[2]) { //2pz2pz
        return pow(1.0 + tau, 2.5) * pow(1.0 - tau, 2.5);
    }
    if (isReversed(a, b)) {
        return calculateOverlapSTOSamePosition(-tau, rho, -kappa, rho_alpha, rho_beta, b, a);
    }
    return 0;
    //throw std::runtime_error("Failed to catch case");
}

double calculateOverlapSTO(double tau, double rho, double kappa, double rho_alpha, double rho_beta, int *a, int *b) {
    if (fabs(tau) < tolerance) {
        return calculateOverlapSTOSameZeta(tau, rho, kappa, rho_alpha, rho_beta, a, b);
    }
    else if (fabs(rho) < tolerance) {
        return calculateOverlapSTOSamePosition(tau, rho, kappa, rho_alpha, rho_beta, a, b);
    }
    else {
        return calculateOverlapSTOFull(tau, rho, kappa, rho_alpha, rho_beta, a, b);
    }
}

double calculateOverlapSTO(STOFunction a, STOFunction b) {
    double r = distance(a.x, a.y, a.z, b.x, b.y, b.z);
    double tau = (a.zeta - b.zeta) / (a.zeta + b.zeta);
    double rho = 0.5 * (a.zeta + b.zeta) * r;
    double kappa = 0.5 * (tau + 1 / tau);
    double rho_alpha = a.zeta * r, rho_beta = b.zeta * r;
    int avec[3] = {a.nlm.n, a.nlm.l, a.nlm.m};
    int bvec[3] = {b.nlm.n, b.nlm.l, b.nlm.m};
    return calculateOverlapSTO(tau, rho, kappa, rho_alpha, rho_beta, avec, bvec);
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

    for (int k = 0; k < 3; k++) {
        for (int i = 0; i < 3; i++) {
            if (isnan(rotationMatrix(k, i))) {
                rotationMatrix = arma::eye(3, 3);
            }
        }
    }
    return rotationMatrix;
}

//[a,1Sb]
double calculateBasicIntegral(double zeta, double rho, int *a) {
    if (a[0] == 1 && a[1] == 0) {
        return (zeta / rho) * (1.0 - (1.0 + rho) * exp(-2.0 * rho));
    }
    else if (a[0] == 2 && a[1] == 0) {
        return (zeta / rho) * (1.0 - (1.0 + 4.0 / 3.0 * rho + 2.0 / 3.0 * pow(rho, 2)) * exp(-2.0 * rho));
    }
    else {
        return 0;
    }
}

//[1Sa,1Sb]
double calculateBasicCoulombIntegralFull(double zeta, double tau, double rho, double kappa, double rho_alpha, double rho_beta, int *a, int *b) {
    if (a[0] == 1 && a[1] == 0 && b[0] == 1 && b[1] == 0) { //1s1s
        return (zeta / rho) * (1.0 - pow(1 - kappa, 2.0) * (1.0 / 4.0 * (2.0 + kappa) + 1.0 / 4.0 * rho_alpha) * exp(-2.0 * rho_alpha)
                               - (1.0 + pow(kappa, 2.0)) * (1.0 / 4.0 * (2.0 - kappa) + 1.0 / 4.0 * rho_beta) * exp(-2.0 * rho_beta));
    }
    if (a[0] == 1 && a[1] == 0 && b[0] == 2 && b[1] == 0) { //1s2s
        return (zeta / rho) * (1.0 - pow(1 - kappa, 2.0) * (1.0 / 4.0 * (1.0 - kappa - pow(kappa, 2.0)) + 1.0 / 12.0 * (1.0 - 2.0 * kappa) * rho_alpha) * exp(-2 * rho_alpha)
                               - (1.0 + pow(kappa, 2.0)) * (1.0 / 4.0 * (3.0 - 3.0 * kappa - pow(kappa, 2.0)) + 1.0 / 3.0 * (2.0 - kappa) * rho_beta + 1.0 / 6.0 * pow(rho_beta, 2)) * exp(-2 * rho_beta));
    }
    if (a[0] == 2 && a[1] == 0 && b[0] == 1 && b[1] == 0) { //1s2s
        return (zeta / rho) * (1.0 - pow(1 - kappa, 2.0) * (1.0 / 4.0 * (1.0 - kappa - pow(kappa, 2.0)) + 1.0 / 12.0 * (1.0 - 2.0 * kappa) * rho_beta) * exp(-2 * rho_beta)
                               - (1.0 + pow(kappa, 2.0)) * (1.0 / 4.0 * (3.0 - 3.0 * kappa - pow(kappa, 2.0)) + 1.0 / 3.0 * (2.0 - kappa) * rho_alpha + 1.0 / 6.0 * pow(rho_alpha, 2)) * exp(-2 * rho_alpha));
    }
    if (a[0] == 2 && a[1] == 0 && b[0] == 2 && b[1] == 0) { //2s2s
        return (zeta / rho) * (1.0 - pow(1 - kappa, 2.0) * (1.0 / 12.0 * (6.0 - kappa - 8 * pow(kappa, 2.0) - 4.0 * pow(kappa, 3.0)) + 1.0 / 3.0 * (1.0 - kappa - pow(kappa, 2)) * rho_alpha + 1.0 / 18.0 * (1 - 2 * kappa) * pow(rho_alpha, 2)) * exp(-2 * rho_alpha)
                               - (1.0 + pow(1 - kappa, 2.0)) * (1.0 / 12.0 * (6.0 + kappa - 8 * pow(kappa, 2.0) + 4.0 * pow(kappa, 3.0)) + 1.0 / 3.0 * (1.0 + kappa + pow(kappa, 2)) * rho_beta + 1.0 / 18.0 * (1 + 2 * kappa) * pow(rho_beta, 2)) * exp(-2 * rho_beta));
    }
    return 0;
}

//[1Sa,1Sb] tau = 0
double calculateBasicCoulombIntegralSameZeta(double zeta, double tau, double rho, double kappa, double rho_alpha, double rho_beta, int *a, int *b) {
    if (a[0] == 1 && a[1] == 0 && b[0] == 1 && b[1] == 0) { //1s1s
        std::cout << " A " << std::endl;
        return (zeta / rho) * (1.0 - (1.0 + 11.0 / 8.0 * rho + 3.0 / 4.0 * pow(rho, 2.0) + 1.0 / 6.0 * pow(rho, 3.0)) * exp(-2.0 * rho));
    }
    if (a[0] == 1 && a[1] == 0 && b[0] == 2 && b[1] == 0) { //1s2s
        std::cout << "b" << std::endl;
        return (zeta / rho) * (1.0 - (1.0 + 71.0 / 48.0 * rho + 23.0 / 24.0 * pow(rho, 2.0) + 1.0 / 3.0 * pow(rho, 3.0) + 1.0 / 18.0 * pow(rho, 3)) * exp(-2.0 * rho));
    }
    if (a[0] == 2 && a[1] == 0 && b[0] == 1 && b[1] == 0) { //1s2s
        std::cout << "b" << std::endl;
        return (zeta / rho) * (1.0 - (1.0 + 71.0 / 48.0 * rho + 23.0 / 24.0 * pow(rho, 2.0) + 1.0 / 3.0 * pow(rho, 3.0) + 1.0 / 18.0 * pow(rho, 3)) * exp(-2.0 * rho));
    }
    if (a[0] == 2 && a[1] == 0 && b[0] == 2 && b[1] == 0) { //2s2s
        std::cout << " c" << std::endl;
        return (zeta / rho) * (1.0 - (1.0 + 37.0 / 24.0 * rho + 13.0 / 12.0 * pow(rho, 2.0) + 4.0 / 9.0 * pow(rho, 3.0) + 1.0 / 9.0 * pow(rho, 3) + 2.0 / 135.0 * pow(rho, 4)) * exp(-2.0 * rho));
    }
    return 0;
}

//[1Sa,1Sb] rho = 0
double calculateBasicCoulombIntegralSamePosition(double zeta, double tau, double rho, double kappa, double rho_alpha, double rho_beta, int *a, int *b) {
    if (a[0] == 1 && a[1] == 0 && b[0] == 1 && b[1] == 0) { //1s1s
        return 1.0 / 8.0 * (1.0 - pow(tau, 2.0)) * (5.0 - pow(tau, 2.0)) * zeta;
    }
    if (a[0] == 1 && a[1] == 0 && b[0] == 2 && b[1] == 0) { //1s2s
        return 1.0 / 48.0 * (1.0 - pow(tau, 2.0)) * (25.0 - 7 * tau - 5 * pow(tau, 2.0) + 3 * pow(tau, 3)) * zeta;
    }
    if (a[0] == 2 && a[1] == 0 && b[0] == 1 && b[1] == 0) { //1s2s
        return 1.0 / 48.0 * (1.0 - pow(-tau, 2.0)) * (25.0 - 7 * -tau - 5 * pow(-tau, 2.0) + 3 * pow(-tau, 3)) * zeta;
    }
    if (a[0] == 2 && a[1] == 0 && b[0] == 2 && b[1] == 0) { //2s2s
        return 1.0 / 24.0 * (1.0 - pow(tau, 2.0)) * (11.0 - 4 * pow(tau, 2.0) + pow(tau, 4)) * zeta;
    }
    return 0;
}

double calculateBasicCoulombIntegral(double zeta, double tau, double rho, double kappa, double rho_alpha, double rho_beta, int *a, int *b) {
    if (std::abs(rho) < tolerance) {
        std::cout << "pos" << std::endl;
        return calculateBasicCoulombIntegralSamePosition(zeta, tau, rho, kappa, rho_alpha, rho_beta, a, b);
    } else if (std::abs(tau) < tolerance) {
        std::cout << "zeta" << std::endl;
        return calculateBasicCoulombIntegralSameZeta(zeta, tau, rho, kappa, rho_alpha, rho_beta, a, b);
    }
    else {
        std::cout << "full" << std::endl;
        return calculateBasicCoulombIntegralFull(zeta, tau, rho, kappa, rho_alpha, rho_beta, a, b);
    }
}


double calculateCoreValenceInteraction(int *a, int *b) { //V_AB

}

double calculateElectrionRepulsionIntegral(int *a, int *b) {

}

} //namespace Semi
