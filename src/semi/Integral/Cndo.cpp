#define ARMA_DONT_USE_WRAPPER
#include "semi/Integral/Cndo.h"

using namespace arma;

namespace Semi {

arma::mat calculateOverlapMatrixGTO(BasisSet<GTOBasis> b) {
    arma::mat Smatrix(b.myBasis.size(), b.myBasis.size());
    for (int k = 0; k < b.myBasis.size(); k++) {
        for (int l = 0; l < b.myBasis.size(); l++) {
            Smatrix(k, l) = Semi::calculateOverlapGTO(b.myBasis[k], b.myBasis[l]);
        }
    }
    return Smatrix;
}


arma::mat calculateOverlapMatrixCGTO(BasisSet<CGTOBasis> b) {
    arma::mat Smatrix(b.myBasis.size(), b.myBasis.size());
    for (int k = 0; k < b.myBasis.size(); k++) {
        for (int l = 0; l < b.myBasis.size(); l++) {
            Smatrix(k, l) = Semi::calculateOverlapCGTO(b.myBasis[k], b.myBasis[l]);
        }
    }
    return Smatrix;
}


arma::mat calculateFockMatrix() {

}

double calculateIonizationPotential(double charge, double l) {
    switch ((int)charge) {
    case 1: return 13.06;
    case 2: return 0;
    case 3: return (l == 0) ? 5.39 : 3.54;
    case 4: return (l == 0) ? 9.32 : 5.96;
    case 5: return (l == 0) ? 14.05 : 8.30;
    case 6: return (l == 0) ? 19.44 : 10.67;
    case 7: return (l == 0) ? 25.58 : 13.19;
    case 8: return (l == 0) ? 32.38 : 15.85;
    case 9: return (l == 0) ? 40.20 : 18.66;
    case 10: return 0;
    default: return 0;
    }
}

std::vector<int> calculateNumberValenceElectrons(double charge) {
    std::vector<int> vals({0, 0});
    if (charge > 4) {
        vals[1] = charge - 4 ;
    }
    if (charge > 2) {
        vals[0] = (charge == 3 ) ? 1 : 2;
    }
    return vals;
}

//U_uu
double calculateCoreHamiltonian(STOBasis a, STOBasis b) {
    double r = distance(a.x, a.y, a.z, b.x, b.y, b.z);
    double zeta_average = 0.5 * (a.zeta + b.zeta);
    double tau = (a.zeta - b.zeta) / (a.zeta + b.zeta);
    double rho = 0.5 * (a.zeta + b.zeta) * r;
    double kappa = 0.5 * (rho + 1 / rho);
    double rho_alpha = a.zeta * r;
    double rho_beta = b.zeta * r;
    int aOrbitalType [3] =  {a.nlm.n, a.nlm.l, a.nlm.m};
    int bOrbitalType [3] = {b.nlm.n, b.nlm.l, b.nlm.m};
    double gamma =  Semi::calculateBasicCoulombIntegral(a.zeta, tau, rho, kappa, rho_alpha, rho_beta, aOrbitalType, bOrbitalType);
    std::vector<int> v = calculateNumberValenceElectrons(a.zeta);
    double U = -(v[0] + v[1] - 1) * gamma - calculateIonizationPotential(a.zeta, a.nlm.l);
    return U;
}

//P_AA
arma::mat calculateTotalChargeDensity() {
    return 0;
}

//P_uv
arma::mat calculateChargeDensity(arma::mat c_v) {
    arma::mat density(c_v.n_cols, c_v.n_cols);
    for (uint u = 0; u < c_v.n_cols; u++) {
        for (uint v = 0; v < c_v.n_cols; v++) {
            for (uint i = 0; i < c_v.n_cols; i++) {
                density(u, v) += 2 * (c_v(i, u) * c_v(i, v));
            }
        }
    }
    return density;
}

//Gamma_uv integral
double calculateElectronRepulsion(STOBasis a, STOBasis b) {
    if (a.nlm.l != 0 || b.nlm.l != 0) {
        return 0;
    }
    else if (a.nlm.n == 1) {
        double r = distance(a.x, a.y, a.z, b.x, b.y, b.z);
        double rho = 0.5 * (a.zeta + a.zeta) * r;
        int aOrbitalType [3] = {a.nlm.n, a.nlm.l, a.nlm.m};
        return a.zeta * Semi::calculateBasicIntegral(a.zeta, rho, aOrbitalType);
    }
    else if (a.nlm.n == 2) {
        double r = distance(a.x, a.y, a.z, b.x, b.y, b.z);
        double rho = 0.5 * (a.zeta + a.zeta) * r;
        int aOrbitalType [3] = {a.nlm.n, a.nlm.l, a.nlm.m};
        return a.zeta * Semi::calculateBasicIntegral(a.zeta, rho, aOrbitalType);
    }
    return 0;
}

//V_AB integral
double calculateNucleurAttraction(STOBasis a, STOBasis b) {
    if (a.nlm.l != 0 || b.nlm.l != 0) {
        return 0;
    }
    else if (a.nlm.n == 1) {
        double r = distance(a.x, a.y, a.z, b.x, b.y, b.z);
        double rho = 0.5 * (a.zeta + a.zeta) * r;
        int aOrbitalType [3] = {a.nlm.n, a.nlm.l, a.nlm.m};
        return a.zeta * Semi::calculateBasicIntegral(a.zeta, rho, aOrbitalType);
    }
    else if (a.nlm.n == 2) {
        double r = distance(a.x, a.y, a.z, b.x, b.y, b.z);
        double rho = 0.5 * (a.zeta + a.zeta) * r;
        int aOrbitalType [3] = {a.nlm.n, a.nlm.l, a.nlm.m};
        return a.zeta * Semi::calculateBasicIntegral(a.zeta, rho, aOrbitalType);
    }
    return 0;
}

//V_AB integral alternative approx
double calculateNucleurAttractionMatagaNishimoto() {
    return 0;
}

//helper method for B_AB
double getBondingParameter(double a) {
    switch ((int)a) {
    case 1: return 9;
    case 2: return 0;
    case 3: return 9;
    case 4: return 13;
    case 5: return 17;
    case 6: return 21;
    case 7: return 25;
    case 8: return 31;
    case 9: return 39;
    case 10: return 0;
    }
    return 0;
}

//B_AB
double calculateBondingParameter(double a, double b) {
    return -0.5 * (getBondingParameter(a) + getBondingParameter(b));
}

//S_uv
arma::mat calculateOverlapMatrix(BasisSet<STOBasis> a) {
    arma::mat Smatrix(a.myBasis.size(), a.myBasis.size());
    for (uint k = 0; k < a.myBasis.size(); k++) {
        for (uint l = 0; l < k; l++) {
            if (k == l) { //diag
                Smatrix(k, l) =  1;
            }
            else if (a.myBasis[k].id == a.myBasis[l].id) { //same atom
                Smatrix(k, l) = 0;
            }
            else {
                double r = distance(a.myBasis[k].x, a.myBasis[k].y, a.myBasis[k].z, a.myBasis[l].x, a.myBasis[l].y, a.myBasis[l].z);
                double zeta_average = 0.5 * (a.myBasis[k].zeta + a.myBasis[l].zeta);
                double tau = (a.myBasis[k].zeta - a.myBasis[l].zeta) / (a.myBasis[k].zeta + a.myBasis[l].zeta);
                double rho = 0.5 * (a.myBasis[k].zeta + a.myBasis[l].zeta) * r;
                double kappa = 0.5 * (rho + 1 / rho);
                double rho_alpha = a.myBasis[k].zeta * r;
                double rho_beta = a.myBasis[l].zeta * r;
                int aOrbitalType [3] =  {a.myBasis[k].nlm.n, a.myBasis[k].nlm.l, a.myBasis[k].nlm.m};
                int bOrbitalType [3] = {a.myBasis[l].nlm.n, a.myBasis[l].nlm.l, a.myBasis[l].nlm.m};
                Smatrix(k, l) = calculateOverlapSTO(tau, rho, kappa, rho_alpha, rho_beta, aOrbitalType, bOrbitalType);
                bool kbool = (k >= 2 && a.myBasis[k - 2].id == a.myBasis[k - 1].id && a.myBasis[k - 1].id == a.myBasis[k].id && a.myBasis[k - 2].nlm.l == a.myBasis[k - 1].nlm.l && a.myBasis[k - 1].nlm.l == a.myBasis[k].nlm.l);
                bool lbool = (l >= 2 && a.myBasis[l - 2].id == a.myBasis[l - 1].id && a.myBasis[l - 1].id == a.myBasis[l].id && a.myBasis[l - 2].nlm.l == a.myBasis[l - 1].nlm.l && a.myBasis[l - 1].nlm.l == a.myBasis[l].nlm.l);
                if (kbool && lbool) {
                    arma::mat rotationMatrix(3, 3);
                    rotationMatrix = findRotation(a.myBasis[k].x, a.myBasis[k].y, a.myBasis[k].z, a.myBasis[l].x, a.myBasis[l].y, a.myBasis[l].z);
                    arma::mat sub(3, 3);
                    sub = Smatrix.submat(k - 2, l - 2, k, l);
                    sub = sub * rotationMatrix;
                    for (int i = 0; i < 3; i++) {
                        for (int j = 0; j < 3; j++) {
                            Smatrix(k - 2 + i, l - 2 + j) = sub(i, j);
                        }
                    }
                }
                else if (kbool) {
                    arma::mat rotationMatrix(3, 3);
                    rotationMatrix = findRotation(a.myBasis[k].x, a.myBasis[k].y, a.myBasis[k].z, a.myBasis[l].x, a.myBasis[l].y, a.myBasis[l].z);
                    arma::mat sub(3, 1);
                    sub = Smatrix.submat(k - 2, l, k, l);
                    sub = arma::trans(arma::trans(sub) * rotationMatrix);
                    for (int i = 0; i < 3; i++) {
                        for (int j = 0; j < 3; j++) {
                            Smatrix(k - 2 + i, l - 2 + j) = sub(i, j);
                        }
                    }
                }
                else if (lbool) {
                    arma::mat rotationMatrix(3, 3);
                    rotationMatrix = findRotation(a.myBasis[k].x, a.myBasis[k].y, a.myBasis[k].z, a.myBasis[l].x, a.myBasis[l].y, a.myBasis[l].z);
                    arma::mat sub(1, 3);
                    sub = Smatrix.submat(k, l - 2, k, l);
                    sub = sub * rotationMatrix;
                    for (int i = 0; i < 3; i++) {
                        for (int j = 0; j < 3; j++) {
                            Smatrix(k - 2 + i, l - 2 + j) = sub(i, j);
                        }
                    }
                }
            }
        }
    }
    for (uint k = 0; k < a.myBasis.size(); k++) {
        for (uint l = k; l < a.myBasis.size(); l++) {
            Smatrix(k, l) = Smatrix(l, k);
        }
    }
    return Smatrix;
}

double distance (double x1, double y1, double z1, double x2, double y2, double z2) {
    return sqrt((pow((x1 - x2), 2) + pow((y1 - y2), 2) + pow((z1 - z2), 2)));
}

} //namespace Semi
