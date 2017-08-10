#define ARMA_DONT_USE_WRAPPER
#include "Cndo.h"
#include "semi/semi_utils.h"
using namespace arma;

namespace Semi {


void SCFold(BasisSet<STOFunction> valenceBasis, arma::mat coefMatrix, arma::mat overlapMatrix, double maxIterations, double convergence, arma::mat &fock) {
    std::vector<arma::mat> occs;
    std::vector<arma::mat> eigvals;

    unsigned int numpi = 0;
    std::vector<int> ids;
    for (int k = 0; k < valenceBasis.myBasis.size(); k++) {
        if (std::find(ids.begin(), ids.end(), valenceBasis.myBasis[k].id) == ids.end()) {
            ids.push_back(valenceBasis.myBasis[k].id);
            numpi += numValence(valenceBasis.myBasis[k].charge);
        }
    }
    unsigned int num_orbitals = numpi / 2;

    arma::mat density;
    (round(100000*overlapMatrix)/100000).print("Overlap Matrix");

    arma::mat previous(coefMatrix.n_rows, coefMatrix.n_rows);
    previous.ones();

    int counter = 0;
    while (counter < maxIterations) {
        //determine occupied coefficients
        arma::mat occ = coefMatrix.cols(0, num_orbitals - 1);
        calculateChargeDensity(occ, density);
        occs.push_back(occ);

        //calculate Fock matrix
        calculateFockMatrix3(valenceBasis, occ, overlapMatrix, fock);

        //solve eigenvalue equation FC = SCE
        arma::mat eigvec;
        arma::vec eigval;
        eig_sym(eigval, eigvec, fock);

        //check convergence
        if (-log10(norm(density - previous)) > convergence) {
            std::cout << "Iteration: " << counter << " Convergence: " << norm(density - previous) << std::endl;
            printEnergy(eigval, num_orbitals, 0);
            break;
        }

        //restart loop using solutions as new coefficients
        coefMatrix = eigvec;
        previous = density;
        counter++;
    }
}

void printEnergy(arma::mat eigval, double num_orbitals, double variant) {
    double energy = 0;
    std::cout << "Orbital Energies" << std::endl;
    for (int i = 0; i < eigval.size(); i++) {
        if (eigval(i) > 0)
            std::cout << "Energy = " << "+" << eigval(i) << std::endl;
        else
            std::cout << "Energy = " << "-" << (-1.00 * eigval(i)) << std::endl;
    }
    for (int i = 0; i < num_orbitals; i++) {
        energy += 2.0 * eigval(i);
    }
    std::cout << "Total Energy: " << energy << std::endl;
}


/********************************************************************************************************************************/
/********************************************************************************************************************************/
/********************************************************************************************************************************/
/********************************************************************************************************************************/
/********************************************************************************************************************************/
/********************************************************************************************************************************/
/********************************************************************************************************************************/
/********************************************************************************************************************************/

void calculateFockMatrix(BasisSet<STOFunction> a, arma::mat occ, arma::mat S, arma::mat &fock) {
    fock.set_size(a.myBasis.size(), a.myBasis.size());
    fock.zeros();

    arma::mat density;
    calculateChargeDensity(occ, density);

    for (int u = 0; u < a.myBasis.size(); u++) {
        for (int v = 0; v < a.myBasis.size(); v++) {
            if (u != v && a.myBasis[u].id != a.myBasis[v].id) {
                fock(u, v) = calculateBondingParameter(a.myBasis[u].charge, a.myBasis[v].charge) * S(u, v)
                             -   0.5 * (density(u, v) * calculateElectronRepulsion(a.myBasis[u], a.myBasis[v]));
            }
            else if (u != v && a.myBasis[u].id == a.myBasis[v].id) {
                fock(u, v) = - 0.5 * (density(u, v) * calculateElectronRepulsion(a.myBasis[u], a.myBasis[u]));
            }
            else {
                fock(u, u) = calculateCoreHamiltonian(a.myBasis[u], a.myBasis[v])
                             + (calculateTotalChargeDensity(density, a.myBasis[u].id, a) - 0.5 * density(u, u)) * calculateElectronRepulsion(a.myBasis[u], a.myBasis[u]);
                std::vector<int> ids;
                double sum = 0;
                for (int k = 0; k < a.myBasis.size(); k++) {
                    if (a.myBasis[k].id != a.myBasis[u].id && std::find(ids.begin(), ids.end(), a.myBasis[k].id) == ids.end()) {
                        ids.push_back(a.myBasis[k].id);
                        fock(u, u) += calculateTotalChargeDensity(density, a.myBasis[k].id, a) * calculateElectronRepulsion(a.myBasis[u], a.myBasis[k])
                                      - calculateNucleurAttraction(a.myBasis[u], a.myBasis[k]);
                    }
                }
            }
        }
    }
}

//U_uu
double calculateCoreHamiltonian(STOFunction a, STOFunction b) {
    double gamma =  calculateElectronRepulsion(a, b);
    std::vector<int> v = calculateNumberValenceElectrons(a.charge);
    return -(v[0] + v[1] - 1) * gamma - calculateIonizationPotential(a.charge, a.nlm.l);
}

std::vector<int> calculateNumberValenceElectrons(double charge) {
    std::vector<int> valence;
    switch (int(charge)) {
    case 1: valence = {1, 0}; break;
    case 2: valence = {2, 0}; break;
    case 3: valence = {1, 0}; break;
    case 4: valence = {2, 0}; break;
    case 5: valence = {2, 1}; break;
    case 6: valence = {2, 2}; break;
    case 7: valence = {2, 3}; break;
    case 8: valence = {2, 4}; break;
    case 9: valence = {2, 5}; break;
    case 10: valence = {2, 6}; break;
    }
    return valence;
}

//U_uu, v_ab
double calculateIonizationPotential(double charge, double l) {
    switch ((int)charge) {
    case 1: return 13.06 / HARTREE_TO_EV;
    case 2: return 0;
    case 3: return (l == 0) ? 5.39 / HARTREE_TO_EV : 3.54 / HARTREE_TO_EV;
    case 4: return (l == 0) ? 9.32 / HARTREE_TO_EV : 5.96 / HARTREE_TO_EV;
    case 5: return (l == 0) ? 14.05 / HARTREE_TO_EV : 8.30 / HARTREE_TO_EV;
    case 6: return (l == 0) ? 19.44 / HARTREE_TO_EV : 10.67 / HARTREE_TO_EV;
    case 7: return (l == 0) ? 25.58 / HARTREE_TO_EV : 13.19 / HARTREE_TO_EV;
    case 8: return (l == 0) ? 32.38 / HARTREE_TO_EV : 15.85 / HARTREE_TO_EV;
    case 9: return (l == 0) ? 40.20 / HARTREE_TO_EV : 18.66 / HARTREE_TO_EV;
    case 10: return 0;
    default: return 0;
    }
}

//P_AA
double calculateTotalChargeDensity(arma::mat density, double id, BasisSet<STOFunction> a) {
    double sum = 0;
    for (int k = 0; k < a.myBasis.size(); k++) {
        if (abs(a.myBasis[k].id - id) < 0.00001) {
            sum += density(k, k);
        }
    }
    return sum;
}

//P_uv
void calculateChargeDensity(arma::mat c_v, arma::mat & density) {
    density.set_size(c_v.n_rows, c_v.n_rows);
    density =  2.0 * c_v * trans(c_v);
}

//Gamma_uv integral
double calculateElectronRepulsion(STOFunction a, STOFunction b) {
    double r = distance(a.x, a.y, a.z, b.x, b.y, b.z);
    double zeta_average = 0.5 * (a.zeta + b.zeta);
    double tau = (a.zeta - b.zeta) / (a.zeta + b.zeta);
    double rho = 0.5 * (a.zeta + b.zeta) * r;
    double kappa = 0.5 * (tau + 1.0 / tau);
    double rho_alpha = a.zeta * r;
    double rho_beta = b.zeta * r;
    double az = 1, bz = 1;

    if (a.charge > 1) {
        az = 2;
    }
    if (b.charge > 1) {
        bz = 2;
    }
    QNumber aQNumber(az, 0, 0);
    QNumber bQNumber(bz, 0, 0);
    double gamma =  Semi::calculateBasicCoulombIntegral(zeta_average, tau, rho, kappa, rho_alpha, rho_beta, aQNumber, bQNumber);
    return gamma;
}

//V_AB integral
double calculateNucleurAttraction(STOFunction a, STOFunction b) {
    double num = b.charge;
    if (b.charge > 2) {
        num -= 2 ;
    }
    double r = distance(a.x, a.y, a.z, b.x, b.y, b.z);
    double rho = a.zeta * r;

    if (a.nlm.n == 1) {
        return (a.zeta / rho) * (1.0 - (1.0 + rho) * exp(-2.0 * rho)) * num;
    }
    else if (a.nlm.n == 2) {
        return (a.zeta / rho) * (1.0 - (1.0 + 3.0 / 2.0 * rho + pow(rho, 2) + 1.0 / 3.0 * pow(rho, 3)) * exp(-2.0 * rho)) * num;
    }
    return 0;
}

//B_AB
double calculateBondingParameter(double a, double b) {
    return -0.5 * (getBondingParameter(a) + getBondingParameter(b));
}

//helper method for B_AB
double getBondingParameter(double a) {
    double h = 27.21138602;
    switch ((int)a) {
    case 1: return 9.0 / HARTREE_TO_EV;
    case 2: return 0;
    case 3: return 9.0 / HARTREE_TO_EV;
    case 4: return 13.0 / HARTREE_TO_EV;
    case 5: return 17.0 / HARTREE_TO_EV;
    case 6: return 21.0 / HARTREE_TO_EV;
    case 7: return 25.0 / HARTREE_TO_EV;
    case 8: return 31.0 / HARTREE_TO_EV;
    case 9: return 39.0 / HARTREE_TO_EV;
    case 10: return 0 / HARTREE_TO_EV;
    }
    return 0;
}

void calculateFockMatrix3(BasisSet<STOFunction> a, arma::mat occ, arma::mat S, arma::mat &fock) {
    fock.set_size(a.myBasis.size(), a.myBasis.size());
    fock.zeros();

    arma::mat density;
    calculateChargeDensity(occ, density);

    for (int u = 0; u < a.myBasis.size(); u++) {
        for (int v = 0; v < a.myBasis.size(); v++) {
            if (u != v && a.myBasis[u].id != a.myBasis[v].id) {
                fock(u, v) = calculateBondingParameter(a.myBasis[u].charge, a.myBasis[v].charge) * S(u, v)
                             -   0.5 * (density(u, v) * calculateElectronRepulsion(a.myBasis[u], a.myBasis[v]));
            }
            else if (u != v && a.myBasis[u].id == a.myBasis[v].id) {
                fock(u, v) = - 0.5 * (density(u, v) * calculateElectronRepulsion(a.myBasis[u], a.myBasis[u]));
            }
            else {
                fock(u, u) = -calculateAtomicData(a.myBasis[u].charge, a.myBasis[u].nlm.l)
                             + (calculateTotalChargeDensity(density, a.myBasis[u].id, a) - numValence(a.myBasis[u].charge) - 0.5 * (density(u, u) - 1)) * calculateElectronRepulsion(a.myBasis[u], a.myBasis[u]);
                std::vector<int> ids;
                for (int k = 0; k < a.myBasis.size(); k++) {
                    if (a.myBasis[k].id != a.myBasis[u].id && std::find(ids.begin(), ids.end(), a.myBasis[k].id) == ids.end()) {
                        ids.push_back(a.myBasis[k].id);
                        fock(u, u) += (calculateTotalChargeDensity(density, a.myBasis[k].id, a) - numValence(a.myBasis[k].charge)) * calculateElectronRepulsion(a.myBasis[u], a.myBasis[k]);
                    }
                }
            }
        }
    }
}

double calculateAtomicData(double charge, double l) {
    switch ((int)charge) {
    case 1: return 7.176 / HARTREE_TO_EV;
    case 2: return 0;
    case 3: return (l == 0) ? 3.106 / HARTREE_TO_EV : 1.258 / HARTREE_TO_EV;
    case 4: return (l == 0) ? 5.946 / HARTREE_TO_EV : 2.563 / HARTREE_TO_EV;
    case 5: return (l == 0) ? 9.594 / HARTREE_TO_EV : 4.001 / HARTREE_TO_EV;
    case 6: return (l == 0) ? 14.051 / HARTREE_TO_EV : 5.572 / HARTREE_TO_EV;
    case 7: return (l == 0) ? 19.316 / HARTREE_TO_EV : 7.275 / HARTREE_TO_EV;
    case 8: return (l == 0) ? 25.390 / HARTREE_TO_EV : 9.111 / HARTREE_TO_EV;
    case 9: return (l == 0) ? 32.272 / HARTREE_TO_EV : 11.080 / HARTREE_TO_EV;
    case 10: return 0;
    default: return 0;
    }
}

//V_AB integral alternative approx
double calculateElectronRepulsionMatagaNishimoto(STOFunction a, STOFunction b) {
    return getIE(a.charge) - getEA(a.charge);
}

double calculateElectronRepulsionMatagaNishimotoOffDiag(STOFunction a, STOFunction b) {
    double x = 2 * pow(2.718, 2) / (calculateElectronRepulsionMatagaNishimoto(a, a) + calculateElectronRepulsionMatagaNishimoto(b, b));
    return pow(2.718, 2) / (distance(a.x, a.y, a.z, b.x, b.y, b.z) + x);
}


















void calculateOverlapMatrixGTO(BasisSet<GTOFunction> b, arma::mat & Smatrix) {
    Smatrix.set_size(b.myBasis.size(), b.myBasis.size());
    for (int k = 0; k < b.myBasis.size(); k++) {
        for (int l = 0; l < b.myBasis.size(); l++) {
            if (k == l) {
                Smatrix(k, l) = 1;
            }
            else {
                Smatrix(k, l) = Semi::calculateOverlapGTO(b.myBasis[k], b.myBasis[l]);
            }
        }
    }
}

void calculateOverlapMatrixSTO(BasisSet<STOFunction> b, arma::mat & Smatrix) {
    Smatrix.set_size(b.myBasis.size(), b.myBasis.size());
    for (int k = 0; k < b.myBasis.size(); k++) {
        for (int l = 0; l < b.myBasis.size(); l++) {
            if (k == l) {
                Smatrix(k, l) = 1;
            }
            else {
                Smatrix(k, l) = Semi::calculateOverlapSTO(b.myBasis[k], b.myBasis[l]);
            }
        }
    }
}

void calculateOverlapMatrixCGTO(BasisSet<CGTOFunction> b, arma::mat & Smatrix) {
    Smatrix.set_size(b.myBasis.size(), b.myBasis.size());
    for (int k = 0; k < b.myBasis.size(); k++) {
        for (int l = 0; l < b.myBasis.size(); l++) {
            if (k == l) {
                Smatrix(k, l) = 1;
            }
            else {
                Smatrix(k, l) = Semi::calculateOverlapCGTO(b.myBasis[k], b.myBasis[l]);
            }
        }
    }
}





} //namespace Semi
