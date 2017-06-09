#define ARMA_DONT_USE_WRAPPER
#include "Cndo.h"

using namespace arma;

namespace Semi {

arma::mat calculateOverlapMatrixGTO(BasisSet<GTOFunction> b) {
    arma::mat Smatrix(b.myBasis.size(), b.myBasis.size());
    for (int k = 0; k < b.myBasis.size(); k++) {
        for (int l = 0; l < b.myBasis.size(); l++) {
            Smatrix(k, l) = Semi::calculateOverlapGTO(b.myBasis[k], b.myBasis[l]);
        }
    }
    return Smatrix;
}

arma::mat calculateOverlapMatrixCGTO(BasisSet<CGTOFunction> b) {
    arma::mat Smatrix(b.myBasis.size(), b.myBasis.size()); Smatrix.fill(0.0);
    for (int k = 0; k < b.myBasis.size(); k++) {
        for (int l = 0; l < b.myBasis.size(); l++) {
            Smatrix(k, l) = Semi::calculateOverlapCGTO(b.myBasis[k], b.myBasis[l]);
        }
    }
    return Smatrix;
}

arma::mat calculateFockMatrix(BasisSet<STOFunction> a, arma::mat coefMatrix) {
    // arma::mat fock(a.myBasis.size(), a.myBasis.size());
    // arma::mat density = calculateChargeDensity(coefMatrix);
    // for (int u = 0; u < a.myBasis.size(); u++) {
    //     for (int v = 0; v < a.myBasis.size(); v++) {
    //         if (u == v) {
    //             fock(u, u) = calculateCoreHamiltonian(a.myBasis[u], a.myBasis[u])
    //                          + (calculateTotalChargeDensity(density, a.myBasis[u].id, a) + 0.5 * density(u, u)) * calculateElectronRepulsion(a.myBasis[u], a.myBasis[u]);
    //             for (int k = 0; k < x; k++) {

    //             }
    //         }
    //         else {
    //             fock(u, v) = calculateBondingParameter(a.myBasis[u].zeta, a.myBasis[v].zeta) * S(u, v)
    //                          - 0.5 * (density(u, v) * calculateElectronRepulsion(a.myBasis[u], a.myBasis[v]));
    //         }
    //     }
    // }

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
double calculateCoreHamiltonian(STOFunction a, STOFunction b) {
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
double calculateTotalChargeDensity(arma::mat density, double id, BasisSet<STOFunction> a) {
    double sum = 0;
    for (int k = 0; k < a.myBasis.size(); k++) {
        if (abs(a.myBasis[k].id == id) < 0.00001) {
            sum += density(k, k);
        }
    }
}

//P_uv
arma::mat calculateChargeDensity(arma::mat c_v) {
    arma::mat density(c_v.n_cols, c_v.n_cols);
    for (unsigned u = 0; u < c_v.n_cols; u++) {
        for (unsigned v = 0; v < c_v.n_cols; v++) {
            for (unsigned i = 0; i < c_v.n_cols; i++) {
                density(u, v) += 2 * (c_v(i, u) * c_v(i, v));
            }
        }
    }
    return density;
}

//Gamma_uv integral
double calculateElectronRepulsion(STOFunction a, STOFunction b) {
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
double calculateNucleurAttraction(STOFunction a, STOFunction b) {
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
arma::mat calculateOverlapMatrix(BasisSet<STOFunction> a) {
    arma::mat Smatrix(a.myBasis.size(), a.myBasis.size());
    for (unsigned k1 = 0; k1 < a.myBasis.size(); k1++) {
        for (unsigned l1 = 0; l1 < a.myBasis.size(); l1++) {
            if (k1 == l1) {
                Smatrix(k1, l1) =  1;
            }
            else if (a.myBasis[k1].id == a.myBasis[l1].id) {
                Smatrix(k1, l1) = calculateOverlapSTO(a.myBasis[k1], a.myBasis[l1]);;
            }
            else {

                Smatrix(k1, l1) = calculateOverlapSTO(a.myBasis[k1], a.myBasis[l1]);
                //Smatrix(k, l) = 3;

            }
        }
    }
    for (unsigned k = 0; k < a.myBasis.size(); k++) {
        for (unsigned l = 0; l < a.myBasis.size(); l++) {

            std::cout << k << l << std::endl;
            std::cout << a.myBasis[k].nlm.l << a.myBasis[l].nlm.l << std::endl;

            if (a.myBasis[k].id == a.myBasis[k - 1].id && a.myBasis[k].nlm.l == a.myBasis[k - 1].nlm.l && a.myBasis[k].id == a.myBasis[k + 1].id && a.myBasis[k].nlm.l == a.myBasis[k + 1].nlm.l && k > 1 && k < a.myBasis.size() - 2
                    && a.myBasis[l].id == a.myBasis[l - 1].id && a.myBasis[l].nlm.l == a.myBasis[l - 1].nlm.l && a.myBasis[l].id == a.myBasis[l + 1].id && a.myBasis[l].nlm.l == a.myBasis[l + 1].nlm.l && l > 1 && l < a.myBasis.size() - 2) {
                std::cout << "kl" << std::endl;
                arma::mat rotationMatrix = findRotation(a.myBasis[k].x, a.myBasis[k].y, a.myBasis[k].z, a.myBasis[l].x, a.myBasis[l].y, a.myBasis[l].z);
                arma::mat sub(3, 3);
                sub = Smatrix.submat(k - 1, l - 1, k + 1, l + 1);
                sub = rotationMatrix * sub  * trans(rotationMatrix);
                Smatrix.print();
                rotationMatrix.print();
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        Smatrix(k - 1 + i, l - 1 + j) = sub(i, j);
                    }
                }
            }
            else if (a.myBasis[l].nlm.l == 0 && a.myBasis[k].id == a.myBasis[k - 1].id && a.myBasis[k].nlm.l == a.myBasis[k - 1].nlm.l && a.myBasis[k].id == a.myBasis[k + 1].id && a.myBasis[k].nlm.l == a.myBasis[k + 1].nlm.l && k > 1 && k < a.myBasis.size() - 2) {
                std::cout << "k" << std::endl;
                arma::mat rotationMatrix = findRotation(a.myBasis[k].x, a.myBasis[k].y, a.myBasis[k].z, a.myBasis[l].x, a.myBasis[l].y, a.myBasis[l].z);
                arma::mat sub(3, 1);
                sub = Smatrix.submat(k - 1, l, k + 1, l);
                sub = rotationMatrix * sub;
                Smatrix.print();
                rotationMatrix.print();
                for (int i = 0; i < 3; i++) {
                    Smatrix(k  + i - 1, l ) = sub(i, 0);
                }
            }
            else if (a.myBasis[k].nlm.l == 0 && a.myBasis[l].id == a.myBasis[l - 1].id && a.myBasis[l].nlm.l == a.myBasis[l - 1].nlm.l && a.myBasis[l].id == a.myBasis[l + 1].id && a.myBasis[l].nlm.l == a.myBasis[l + 1].nlm.l && l > 1 && l < a.myBasis.size() - 2) {
                std::cout << "l" << std::endl;
                arma::mat rotationMatrix = findRotation(a.myBasis[k].x, a.myBasis[k].y, a.myBasis[k].z, a.myBasis[l].x, a.myBasis[l].y, a.myBasis[l].z);
                arma::mat sub(1, 3);
                sub = trans(Smatrix.submat(k, l - 1, k, l + 1));
                sub = trans(rotationMatrix * sub);
                Smatrix.print();
                rotationMatrix.print();
                for (int i = 0; i < 3; i++) {
                    Smatrix(k, l   + i - 1) = sub(0, i);
                }
            }
        }
    }
    return Smatrix;
}

} //namespace Semi
