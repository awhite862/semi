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

arma::mat calculateOverlapMatrixSTO(BasisSet<STOFunction> b) {
    arma::mat Smatrix(b.myBasis.size(), b.myBasis.size());
    for (int k = 0; k < b.myBasis.size(); k++) {
        for (int l = 0; l < b.myBasis.size(); l++) {
            Smatrix(k, l) = Semi::calculateOverlapSTO(b.myBasis[k], b.myBasis[l]);
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

void SCF(BasisSet<STOFunction> a, arma::mat coefMatrix, arma::mat S, arma::mat &fock) {
    std::vector<arma::mat> temp;
    std::vector<arma::mat> occs;
    std::vector<arma::mat> f;
    for (int k = 0; k < 5; k++) {
        calculateFockMatrix(a, coefMatrix, S, fock);
        arma::mat occ = coefMatrix.cols(0, 3);
        arma::mat density = calculateChargeDensity(occ);
        f.push_back(fock);
        occs.push_back(occ);
        temp.push_back(coefMatrix);
        arma::mat eigvec;
        arma::vec eigval;
        eig_sym(eigval, eigvec, fock);
        coefMatrix = eigvec;
    }
    (temp[0]).print();
    std::cout << "fock" << std::endl;
    f[0].print();
    for (int k = 1; k < 5; k++) {
        std::cout << k << " " << norm(calculateChargeDensity(temp[k]) - calculateChargeDensity(temp[k - 1])) << std::endl;
        std::cout << k << " " << norm(calculateChargeDensity(occs[k]) - calculateChargeDensity(occs[k - 1])) << std::endl;
        (temp[k]).print();
        std::cout << "fock" << std::endl;
        f[k].print();
        std::cout << "trace" << trace(calculateChargeDensity(occs[k])) << std::endl;
    }
}

//tr(density) = num_elec

void calculateFockMatrix(BasisSet<STOFunction> a, arma::mat coefMatrix, arma::mat S, arma::mat &fock) {
    fock.zeros();
    arma::mat occ = coefMatrix.cols(0, 3);
    arma::mat density = calculateChargeDensity(occ);
    coefMatrix.print("coeffs");
    occ.print("occ");
    density.print("density");
    for (int u = 0; u < a.myBasis.size(); u++) {
        for (int v = 0; v < a.myBasis.size(); v++) {
            fock.print("f");
            if (u == v) {
                std::cout << "diag fock:" << fock(u, u) << std::endl;

                std::cout << u << v << " eq" << std::endl;
                std::cout << "core: " << calculateCoreHamiltonian(a.myBasis[u], a.myBasis[u]) << std::endl;
                std::cout << "charge density: " << calculateTotalChargeDensity(density, a.myBasis[u].id, a) << std::endl;
                std::cout << "density: " << density(u, u) << std::endl;
                std::cout << "electron repulsion: " << calculateElectronRepulsion(a.myBasis[u], a.myBasis[u]) << std::endl;
                fock(u, u) = calculateCoreHamiltonian(a.myBasis[u], a.myBasis[u])
                             + (calculateTotalChargeDensity(density, a.myBasis[u].id, a) + 0.5 * density(u, u)) * calculateElectronRepulsion(a.myBasis[u], a.myBasis[u]);
                std::vector<int> ids;
                std::cout << "diag fock:" << fock(u, u) << std::endl;
                for (int k = 0; k < a.myBasis.size(); k++) {
                    if (a.myBasis[k].id != a.myBasis[u].id && std::find(ids.begin(), ids.end(), a.myBasis[k].id) == ids.end()) {
                        ids.push_back(a.myBasis[k].id);
                        std::cout << "loop" << std::endl;
                        std::cout << "charge density: " << calculateTotalChargeDensity(density, a.myBasis[k].id, a) << std::endl;
                        std::cout << "electron repulsion: " << calculateElectronRepulsion(a.myBasis[u], a.myBasis[k]) << std::endl;
                        std::cout << "nucl attraction: " << calculateNucleurAttraction(a.myBasis[u], a.myBasis[k]) << std::endl;
                        fock(u, u) += calculateTotalChargeDensity(density, a.myBasis[k].id, a) * calculateElectronRepulsion(a.myBasis[u], a.myBasis[k])
                                      - calculateNucleurAttraction(a.myBasis[u], a.myBasis[k]);
                        std::cout << "diag fock:" << fock(u, u) << std::endl;
                    }
                }
            }
            else {
                std::cout << u << v << " neq" << std::endl;
                std::cout << "bonding:" << calculateBondingParameter(a.myBasis[u].zeta, a.myBasis[v].zeta) << " overlap: " <<  S(u, v) << " density: " << density(u, v) << " e- repuls: " << calculateElectronRepulsion(a.myBasis[u], a.myBasis[v]) << std::endl;
                fock(u, v) = calculateBondingParameter(a.myBasis[u].zeta, a.myBasis[v].zeta) * S(u, v)
                             - 0.5 * (density(u, v) * calculateElectronRepulsion(a.myBasis[u], a.myBasis[v]));
            }
        }
    }
}


//diag should be negativ
//U_uu, v_ab
double calculateIonizationPotential(double charge, double l) {
    double h = 27.21138602;
    switch ((int)charge) {
    case 1: return 13.06 / h;
    case 2: return 0;
    case 3: return (l == 0) ? 5.39 / h : 3.54 / h;
    case 4: return (l == 0) ? 9.32 / h : 5.96 / h;
    case 5: return (l == 0) ? 14.05 / h : 8.30 / h;
    case 6: return (l == 0) ? 19.44 / h : 10.67 / h;
    case 7: return (l == 0) ? 25.58 / h : 13.19 / h;
    case 8: return (l == 0) ? 32.38 / h : 15.85 / h;
    case 9: return (l == 0) ? 40.20 / h : 18.66 / h;
    case 10: return 0;
    default: return 0;
    }
}

std::vector<int> calculateNumberValenceElectrons(double charge) {
    std::vector<int> vals({1, 0});
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
    double gamma =  Semi::calculateBasicCoulombIntegral(a.zeta, b.zeta, tau, rho, kappa, rho_alpha, rho_beta, aOrbitalType, bOrbitalType);
    std::vector<int> v = calculateNumberValenceElectrons(a.zeta);
    //std::cout << "valence" << a.zeta << " " << v[0] << " " << v[1] << std::endl;
    double U = -(v[0] + v[1] - 1) * gamma - calculateIonizationPotential(a.zeta, a.nlm.l);
    return U;
}

//P_AA
double calculateTotalChargeDensity(arma::mat density, double id, BasisSet<STOFunction> a) {
    double sum = 0;
    //std::vector<int> ids;
    for (int k = 0; k < a.myBasis.size(); k++) {
        //if (std::find(ids.begin(), ids.end(), a.myBasis[k].id) == ids.end()) {
        if (abs(a.myBasis[k].id - id) < 0.00001) {
            std::cout << k << std::endl;
            //ids.push_back(a.myBasis[k].id);
            sum += density(k, k);
        }
    }

    return sum;
}

//P_uv
arma::mat calculateChargeDensity(arma::mat c_v) {
    // arma::mat density(c_v.n_cols, c_v.n_cols);
    // for (unsigned u = 0; u < c_v.n_cols; u++) {
    //     for (unsigned v = 0; v < c_v.n_cols; v++) {
    //         for (unsigned i = 0; i < c_v.n_cols; i++) {
    //             density(u, v) += 2 * (c_v(i, u) * c_v(i, v));
    //         }
    //     }
    // }
    // return density;

    return 2 * c_v * trans(c_v);
}

//Gamma_uv integral
double calculateElectronRepulsion(STOFunction a, STOFunction b) {
    double r = distance(a.x, a.y, a.z, b.x, b.y, b.z);
    double zeta_average = 0.5 * (a.zeta + b.zeta);
    double tau = (a.zeta - b.zeta) / (a.zeta + b.zeta);
    double rho = 0.5 * (a.zeta + b.zeta) * r;
    double kappa = 0.5 * (rho + 1 / rho);
    double rho_alpha = a.zeta * r;
    double rho_beta = b.zeta * r;
    int aOrbitalType [3] =  {a.nlm.n, a.nlm.l, a.nlm.m};
    int bOrbitalType [3] = {b.nlm.n, b.nlm.l, b.nlm.m};
    double gamma =  Semi::calculateBasicCoulombIntegral(a.zeta, b.zeta, tau, rho, kappa, rho_alpha, rho_beta, aOrbitalType, bOrbitalType);
    return gamma;
}

//V_AB integral
double calculateNucleurAttraction(STOFunction a, STOFunction b) {
    double r = distance(a.x, a.y, a.z, b.x, b.y, b.z);
    double rho = a.zeta * r;
    int aOrbitalType [3] = {a.nlm.n, a.nlm.l, a.nlm.m};
    if (aOrbitalType[0] == 1 && aOrbitalType[1] == 0) {
        return (a.zeta / rho) * (1.0 - (1.0 + rho) * exp(-2.0 * rho));
    }
    else if (aOrbitalType[0] == 2 && aOrbitalType[1] == 0) {
        return (a.zeta / rho) * (1.0 - (1.0 + 4.0 / 3.0 * rho + 2.0 / 3.0 * pow(rho, 2)) * exp(-2.0 * rho));
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
                Smatrix(k1, l1) = calculateOverlapSTO(a.myBasis[k1], a.myBasis[l1]);
            }
            else {
                Smatrix(k1, l1) = calculateOverlapSTO(a.myBasis[k1], a.myBasis[l1]);
            }
        }
    }
    for (unsigned k = 0; k < a.myBasis.size(); k++) {
        for (unsigned l = 0; l < a.myBasis.size(); l++) {
            if (a.myBasis[k].id == a.myBasis[k - 1].id && a.myBasis[k].nlm.l == a.myBasis[k - 1].nlm.l && a.myBasis[k].id == a.myBasis[k + 1].id && a.myBasis[k].nlm.l == a.myBasis[k + 1].nlm.l && k > 1 && k < a.myBasis.size() - 2
                    && a.myBasis[l].id == a.myBasis[l - 1].id && a.myBasis[l].nlm.l == a.myBasis[l - 1].nlm.l && a.myBasis[l].id == a.myBasis[l + 1].id && a.myBasis[l].nlm.l == a.myBasis[l + 1].nlm.l && l > 1 && l < a.myBasis.size() - 2) {
                arma::mat rotationMatrix = findRotation(a.myBasis[k].x, a.myBasis[k].y, a.myBasis[k].z, a.myBasis[l].x, a.myBasis[l].y, a.myBasis[l].z);
                arma::mat sub(3, 3);
                sub = Smatrix.submat(k - 1, l - 1, k + 1, l + 1);
                sub = rotationMatrix * sub  * trans(rotationMatrix);
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        Smatrix(k - 1 + i, l - 1 + j) = sub(i, j);
                    }
                }
            }
            else if (a.myBasis[l].nlm.l == 0 && a.myBasis[k].id == a.myBasis[k - 1].id && a.myBasis[k].nlm.l == a.myBasis[k - 1].nlm.l && a.myBasis[k].id == a.myBasis[k + 1].id && a.myBasis[k].nlm.l == a.myBasis[k + 1].nlm.l && k > 1 && k < a.myBasis.size() - 2) {
                arma::mat rotationMatrix = findRotation(a.myBasis[k].x, a.myBasis[k].y, a.myBasis[k].z, a.myBasis[l].x, a.myBasis[l].y, a.myBasis[l].z);
                arma::mat sub(3, 1);
                sub = Smatrix.submat(k - 1, l, k + 1, l);
                sub = rotationMatrix * sub;
                for (int i = 0; i < 3; i++) {
                    Smatrix(k  + i - 1, l ) = sub(i, 0);
                }
            }
            else if (a.myBasis[k].nlm.l == 0 && a.myBasis[l].id == a.myBasis[l - 1].id && a.myBasis[l].nlm.l == a.myBasis[l - 1].nlm.l && a.myBasis[l].id == a.myBasis[l + 1].id && a.myBasis[l].nlm.l == a.myBasis[l + 1].nlm.l && l > 1 && l < a.myBasis.size() - 2) {
                arma::mat rotationMatrix = findRotation(a.myBasis[k].x, a.myBasis[k].y, a.myBasis[k].z, a.myBasis[l].x, a.myBasis[l].y, a.myBasis[l].z);
                arma::mat sub(1, 3);
                sub = trans(Smatrix.submat(k, l - 1, k, l + 1));
                sub = trans(rotationMatrix * sub);
                for (int i = 0; i < 3; i++) {
                    Smatrix(k, l   + i - 1) = sub(0, i);
                }
            }
        }
    }
    return Smatrix;
}

} //namespace Semi
