#define ARMA_DONT_USE_WRAPPER
#include "Cndo.h"

using namespace arma;

namespace Semi {

void calculateOverlapMatrixGTO(BasisSet<GTOFunction> b, arma::mat &Smatrix) {
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

void calculateOverlapMatrixSTO(BasisSet<STOFunction> b, arma::mat &Smatrix) {
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

void calculateOverlapMatrixCGTO(BasisSet<CGTOFunction> b, arma::mat &Smatrix) {
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

void SCF(BasisSet<STOFunction> a, arma::mat coefMatrix, arma::mat S, arma::mat &fock) {
    std::vector<arma::mat> occs;
    std::vector<arma::mat> f;
    std::vector<arma::mat> eigvals;


    unsigned int numpi = 0;
    std::vector<int> ids;
    for (int k = 0; k < a.myBasis.size(); k++) {
        if (std::find(ids.begin(), ids.end(), a.myBasis[k].id) == ids.end()) {
            ids.push_back(a.myBasis[k].id);
            std::cout << "valence elec" <<  numValence(a.myBasis[k].charge) << std::endl;
            numpi += numValence(a.myBasis[k].charge);
        }
    }
    unsigned int num_orbitals = numpi / 2;

    std::cout << num_orbitals << " orbitals" << std::endl;
    arma::mat density;
    S.print("overlap");
    for (int k = 0; k < 50; k++) {
        calculateFockMatrix3(a, coefMatrix, S, fock);
        arma::mat occ = coefMatrix.cols(0, num_orbitals - 1);
        calculateChargeDensity(occ, density);
        f.push_back(fock);
        occs.push_back(occ);
        arma::mat eigvec;
        arma::vec eigval;
        eig_sym(eigval, eigvec, fock);
        fock.print("fock");
        eigvec.print("eigvecs");

        eigval.print("eigvals");
        eigvals.push_back(eigval);
        coefMatrix = eigvec;


    }
    std::cout << "num orbitals:" << num_orbitals << std::endl;
    std::cout << "-----------------------------------------------------------------------" << std::endl;
    std::cout << "fock" << std::endl;
    f[0].print();
    std::cout << "density occ" << std::endl;
    calculateChargeDensity(occs[0], density);
    density.print();
    arma::mat old;
    double energy = 0;
    for (int k = 1; k < 50; k++) {
        energy = 0;
        std::cout << "-----------------------------------------------------------------------" << std::endl;
        std::cout << "fock" << std::endl;
        (round(1000 * f[k]) / 1000).print();
        old = density;
        std::cout << "occs charge density" << std::endl;
        calculateChargeDensity(occs[k], density);
        (round(1000 * density) / 1000).print();
        std::cout << "iteration: " << k << " convergence: " << norm(old - density) << std::endl;





        std::cout << "Orbital Energies" << std::endl; arma::mat temp = eigvals[k];
        for (int i = 0; i < eigvals[k].size(); i++) {
            temp = eigvals[k];
            if (temp(i) > 0)
                std::cout << "Energy = " << "+" << temp(i) << std::endl;
            else
                std::cout << "Energy = " << "-" << (-1.00 * temp(i)) << std::endl;
        }
        for (int i = 0; i < num_orbitals; i++) {
            temp = eigvals[k];
            energy += 2.0 * temp(i);
        }
        temp.print("eigvals");
        std::cout << "Total Energy: " << energy << std::endl;
    }

}

//tr(density) = num_elec

void calculateFockMatrix(BasisSet<STOFunction> a, arma::mat coefMatrix, arma::mat S, arma::mat & fock) {
    fock.set_size(a.myBasis.size(), a.myBasis.size());
    fock.zeros();
    unsigned int numpi = 0;
    std::vector<int> ids;
    for (int k = 0; k < a.myBasis.size(); k++) {
        if (std::find(ids.begin(), ids.end(), a.myBasis[k].id) == ids.end()) {
            ids.push_back(a.myBasis[k].id);
            std::cout << numValence(a.myBasis[k].charge) << std::endl;
            numpi += numValence(a.myBasis[k].charge);
        }
    }
    unsigned int num_orbitals = numpi / 2;
    std::cout << num_orbitals << std::endl;
    coefMatrix.print();
    arma::mat occ = coefMatrix.cols(0, num_orbitals - 1);
    arma::mat density;
    calculateChargeDensity(occ, density);
    coefMatrix.print("coeffs");
    occ.print("occ");
    density.print("density");

    //debug
    arma::mat bondingParameter(a.myBasis.size(), a.myBasis.size());
    arma::mat electronRepulsion(a.myBasis.size(), a.myBasis.size());
    arma::mat coreHamiltonian(a.myBasis.size(), a.myBasis.size());
    arma::mat chargeDensity(a.myBasis.size(), a.myBasis.size());
    arma::mat second(a.myBasis.size(), a.myBasis.size());
    arma::mat third(a.myBasis.size(), a.myBasis.size());
    arma::mat nucleurAttraction(a.myBasis.size(), a.myBasis.size());
    arma::mat summation(a.myBasis.size(), a.myBasis.size());
    //debug

    bondingParameter.zeros();
    electronRepulsion.zeros();
    coreHamiltonian.zeros();
    chargeDensity.zeros();
    second.zeros();
    third.zeros();
    summation.zeros();
    nucleurAttraction.zeros();

    std::cout << "debug" << std::endl;
    for (int u = 0; u < a.myBasis.size(); u++) {
        for (int v = 0; v < a.myBasis.size(); v++) {
            if (u != v && a.myBasis[u].id != a.myBasis[v].id) {
                fock(u, v) = calculateBondingParameter(a.myBasis[u].charge, a.myBasis[v].charge) * S(u, v)
                             -   0.5 * (density(u, v) * 0.4342/*calculateElectronRepulsion(a.myBasis[u], a.myBasis[v])*/);
                bondingParameter(u, v) = calculateBondingParameter(a.myBasis[u].charge, a.myBasis[v].charge);
                electronRepulsion(u, v) = calculateElectronRepulsion(a.myBasis[u], a.myBasis[v]);
            }
            else if (u != v && a.myBasis[u].id == a.myBasis[v].id) {
                fock(u, v) = - 0.5 * (density(u, v) * calculateElectronRepulsion(a.myBasis[u], a.myBasis[u]));
                electronRepulsion(u, v) = calculateElectronRepulsion(a.myBasis[u], a.myBasis[u]);
            }
            else {
                coreHamiltonian(u, v) = calculateCoreHamiltonian(a.myBasis[u], a.myBasis[u]);
                chargeDensity(u, v) = calculateTotalChargeDensity(density, a.myBasis[u].id, a);
                second(u, v) =  (calculateTotalChargeDensity(density, a.myBasis[u].id, a) - 0.5 * density(u, u)) * calculateElectronRepulsion(a.myBasis[u], a.myBasis[u]);
                fock(u, u) = calculateCoreHamiltonian(a.myBasis[u], a.myBasis[v])
                             + (calculateTotalChargeDensity(density, a.myBasis[u].id, a) - 0.5 * density(u, u)) * calculateElectronRepulsion(a.myBasis[u], a.myBasis[u]);
                std::vector<int> ids;
                electronRepulsion(u, v) = calculateElectronRepulsion(a.myBasis[u], a.myBasis[v]);
                double sum = 0;
                for (int k = 0; k < a.myBasis.size(); k++) {
                    if (a.myBasis[k].id != a.myBasis[u].id && std::find(ids.begin(), ids.end(), a.myBasis[k].id) == ids.end()) {
                        ids.push_back(a.myBasis[k].id);
                        fock(u, u) += calculateTotalChargeDensity(density, a.myBasis[k].id, a) * calculateElectronRepulsion(a.myBasis[u], a.myBasis[k])
                                      - calculateNucleurAttraction(a.myBasis[u], a.myBasis[k]);




                        sum += calculateTotalChargeDensity(density, a.myBasis[k].id, a) * calculateElectronRepulsion(a.myBasis[u], a.myBasis[k])
                               - calculateNucleurAttraction(a.myBasis[u], a.myBasis[k]);
                        summation(u, u) = calculateTotalChargeDensity(density, a.myBasis[k].id, a) * calculateElectronRepulsion(a.myBasis[u], a.myBasis[k]);
                        nucleurAttraction(u, k) = calculateNucleurAttraction(a.myBasis[u], a.myBasis[k]);
                    }
                }
                third(u, v) = sum;
            }
        }
    }
    bondingParameter.print("bonding param");
    electronRepulsion.print("elec repulsion");
    coreHamiltonian.print("core");
    chargeDensity.print("charge density");
    second.print("second term");
    third.print("third term");
    summation.print("summation");
    nucleurAttraction.print("nuc attraction");
    density.print("density");
    S.print("overlap");
}

void calculateFockMatrix3(BasisSet<STOFunction> a, arma::mat coefMatrix, arma::mat S, arma::mat & fock) {
    fock.set_size(a.myBasis.size(), a.myBasis.size());
    fock.zeros();
    unsigned int numpi = 0;
    std::vector<int> ids;
    for (int k = 0; k < a.myBasis.size(); k++) {
        if (std::find(ids.begin(), ids.end(), a.myBasis[k].id) == ids.end()) {
            ids.push_back(a.myBasis[k].id);
            std::cout << numValence(a.myBasis[k].charge) << std::endl;
            numpi += numValence(a.myBasis[k].charge);
        }
    }
    unsigned int num_orbitals = numpi / 2;
    std::cout << num_orbitals << std::endl;
    coefMatrix.print();
    arma::mat occ = coefMatrix.cols(0, num_orbitals - 1);
    arma::mat density;
    calculateChargeDensity(occ, density);
    coefMatrix.print("coeffs");
    occ.print("occ");
    density.print("density");

    //debug
    arma::mat bondingParameter(a.myBasis.size(), a.myBasis.size());
    arma::mat electronRepulsion(a.myBasis.size(), a.myBasis.size());
    arma::mat atomicData(a.myBasis.size(), a.myBasis.size());
    arma::mat second(a.myBasis.size(), a.myBasis.size());
    arma::mat third(a.myBasis.size(), a.myBasis.size());
    arma::mat chargeDensity(a.myBasis.size(), a.myBasis.size());
    arma::mat nucleurAttraction(a.myBasis.size(), a.myBasis.size());
    arma::mat valence(a.myBasis.size(), a.myBasis.size());
    //debug

    bondingParameter.zeros();
    electronRepulsion.zeros();
    atomicData.zeros();
    second.zeros();
    third.zeros();
    chargeDensity.zeros();
    nucleurAttraction.zeros();
    valence.zeros();
    std::cout << "debug" << std::endl;
    for (int u = 0; u < a.myBasis.size(); u++) {
        for (int v = 0; v < a.myBasis.size(); v++) {
            if (u != v && a.myBasis[u].id != a.myBasis[v].id) {
                fock(u, v) = calculateBondingParameter(a.myBasis[u].charge, a.myBasis[v].charge) * S(u, v)
                             -   0.5 * (density(u, v) * calculateElectronRepulsion(a.myBasis[u], a.myBasis[v]));
                bondingParameter(u, v) = calculateBondingParameter(a.myBasis[u].charge, a.myBasis[v].charge);
                electronRepulsion(u, v) = calculateElectronRepulsion(a.myBasis[u], a.myBasis[v]);
            }
            else if (u != v && a.myBasis[u].id == a.myBasis[v].id) {
                fock(u, v) = - 0.5 * (density(u, v) * calculateElectronRepulsion(a.myBasis[u], a.myBasis[u]));
                electronRepulsion(u, v) = calculateElectronRepulsion(a.myBasis[u], a.myBasis[u]);
            }
            else {
                fock(u, u) = -calculateAtomicData(a.myBasis[u].charge, a.myBasis[u].nlm.l)
                             + (calculateTotalChargeDensity(density, a.myBasis[u].id, a) - numValence(a.myBasis[u].charge) - 0.5 * (density(u, u) - 1)) * calculateElectronRepulsion(a.myBasis[u], a.myBasis[u]);
                std::vector<int> ids;
                std::cout << "numvalence:" << u << " " << a.myBasis[u].charge << " " << numValence(a.myBasis[u].charge) << std::endl;
                valence(u, u) = numValence(a.myBasis[u].charge);

                atomicData(u, v) = calculateAtomicData(a.myBasis[u].charge, a.myBasis[u].nlm.l);
                second(u, v) = (calculateTotalChargeDensity(density, a.myBasis[u].id, a) - numValence(a.myBasis[u].charge) - 0.5 * (density(u, u) - 1)) * calculateElectronRepulsion(a.myBasis[u], a.myBasis[u]);
                electronRepulsion(u, v) = calculateElectronRepulsion(a.myBasis[u], a.myBasis[v]);
                for (int k = 0; k < a.myBasis.size(); k++) {
                    if (a.myBasis[k].id != a.myBasis[u].id && std::find(ids.begin(), ids.end(), a.myBasis[k].id) == ids.end()) {
                        ids.push_back(a.myBasis[k].id);
                        third(u, u) = (calculateTotalChargeDensity(density, a.myBasis[k].id, a) - numValence(a.myBasis[k].charge)) * calculateElectronRepulsion(a.myBasis[u], a.myBasis[k]);
                        fock(u, u) += (calculateTotalChargeDensity(density, a.myBasis[k].id, a) - numValence(a.myBasis[k].charge)) * calculateElectronRepulsion(a.myBasis[u], a.myBasis[k]);
                        nucleurAttraction(u, k) = calculateElectronRepulsion(a.myBasis[u], a.myBasis[k]);
                    }
                }
            }
        }
    }
    bondingParameter.print("bonding param");
    electronRepulsion.print("elec repulsion");
    atomicData.print("atomicData");
    chargeDensity.print("charge density");
    second.print("second");
    third.print("third");
    nucleurAttraction.print("gamma");
    valence.print("num valence");
    density.print("density");
    S.print("overlap");
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
    std::vector<int> vals;
    switch (int(charge)) {
    case 1: vals = {1, 0}; break;
    case 2: vals = {2, 0}; break;
    case 3: vals = {1, 0}; break;
    case 4: vals = {2, 0}; break;
    case 5: vals = {2, 1}; break;
    case 6: vals = {2, 2}; break;
    case 7: vals = {2, 3}; break;
    case 8: vals = {2, 4}; break;
    case 9: vals = {2, 5}; break;
    case 10: vals = {2, 6}; break;
    }
    return vals;
}

//U_uu
double calculateCoreHamiltonian(STOFunction a, STOFunction b) {
    double r = distance(a.x, a.y, a.z, b.x, b.y, b.z);
    double zeta_average = 0.5 * (a.zeta + b.zeta);
    double tau = (a.zeta - b.zeta) / (a.zeta + b.zeta);
    double rho = 0.5 * (a.zeta + b.zeta) * r;
    double kappa = 0.5 * (tau + 1 / tau);
    double rho_alpha = a.zeta * r;
    double rho_beta = b.zeta * r;
    int aOrbitalType [3] =  {a.nlm.n, a.nlm.l, a.nlm.m};
    int bOrbitalType [3] = {b.nlm.n, b.nlm.l, b.nlm.m};

    double gamma =  calculateElectronRepulsion(a, b);
    std::vector<int> v = calculateNumberValenceElectrons(a.charge);
    double U = -(v[0] + v[1] - 1) * gamma - calculateIonizationPotential(a.charge, a.nlm.l);
    return U;
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

double calculateAtomicData(double charge, double l) {
    double h = 27.21138602;
    switch ((int)charge) {
    case 1: return 7.176 / h;
    case 2: return 0;
    case 3: return (l == 0) ? 3.106 / h : 1.258 / h;
    case 4: return (l == 0) ? 5.946 / h : 2.563 / h;
    case 5: return (l == 0) ? 9.594 / h : 4.001 / h;
    case 6: return (l == 0) ? 14.051 / h : 5.572 / h;
    case 7: return (l == 0) ? 19.316 / h : 7.275 / h;
    case 8: return (l == 0) ? 25.390 / h : 9.111 / h;
    case 9: return (l == 0) ? 32.272 / h : 11.080 / h;
    case 10: return 0;
    default: return 0;
    }
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
    int aOrbitalType [3] =  {az, 0, 0};
    int bOrbitalType [3] = {bz, 0, 0};
    double gamma =  Semi::calculateBasicCoulombIntegral(zeta_average, zeta_average, tau, rho, kappa, rho_alpha, rho_beta, aOrbitalType, bOrbitalType);
    return gamma;
}

//V_AB integral alternative approx
double calculateElectronRepulsionMatagaNishimoto(STOFunction a, STOFunction b) {
    return getIE(a.charge) - getEA(a.charge);
}

double calculateElectronRepulsionMatagaNishimotoOffDiag(STOFunction a, STOFunction b) {
    double x = 2 * pow(2.718, 2) / (calculateElectronRepulsionMatagaNishimoto(a, a) + calculateElectronRepulsionMatagaNishimoto(b, b));
    return pow(2.718, 2) / (distance(a.x, a.y, a.z, b.x, b.y, b.z) + x);
}

//V_AB integral
double calculateNucleurAttraction(STOFunction a, STOFunction b) {
    double num = b.charge;
    if (b.charge > 2) {
        num -= 2 ;
    }
    double r = distance(a.x, a.y, a.z, b.x, b.y, b.z);
    double rho = a.zeta * r;
    int aOrbitalType [3] = {a.nlm.n, a.nlm.l, a.nlm.m};
    if (aOrbitalType[0] == 1) {
        return (a.zeta / rho) * (1.0 - (1.0 + rho) * exp(-2.0 * rho)) * num;
    }
    else if (aOrbitalType[0] == 2) {
        return (a.zeta / rho) * (1.0 - (1.0 + 3.0 / 2.0 * rho + pow(rho, 2) + 1.0 / 3.0 * pow(rho, 3)) * exp(-2.0 * rho)) * num;
    }
    return 0;
}

//helper method for calculateElectronRepulsionMatagaNishimoto
double getIE(double a) {
    double h = 27.21138602;
    switch ((int)a) {
    case 1: return 13.5984 / h;
    case 2: return 24.5874 / h;
    case 3: return 5.3917 / h;
    case 4: return 9.3227 / h;
    case 5: return 8.298 / h;
    case 6: return 11.2603 / h;
    case 7: return 14.5341 / h;
    case 8: return 13.6181 / h;
    case 9: return 17.4228 / h;
    case 10: return 0;
    }
    return 0;
}

//helper method for calculateElectronRepulsionMatagaNishimoto
double getEA(double a) {
    double h = 27.21138602;
    switch ((int)a) {
    case 1: return 0.754 / h;
    case 2: return -19.7 / h;
    case 3: return 0.618 / h;
    case 4: return -2.4 / h;
    case 5: return 0.279 / h;
    case 6: return 1.262 / h;
    case 7: return -1.4 / h;
    case 8: return 1.461 / h;
    case 9: return 3.401 / h;
    case 10: return 0;
    }
    return 0;
}

//helper method for B_AB
double getBondingParameter(double a) {
    double h = 27.21138602;
    switch ((int)a) {
    case 1: return 9.0 / h;
    case 2: return 0;
    case 3: return 9.0 / h;
    case 4: return 13.0 / h;
    case 5: return 17.0 / h;
    case 6: return 21.0 / h;
    case 7: return 25.0 / h;
    case 8: return 31.0 / h;
    case 9: return 39.0 / h;
    case 10: return 0 / h;
    }
    return 0;
}

//B_AB
double calculateBondingParameter(double a, double b) {
    return -0.5 * (getBondingParameter(a) + getBondingParameter(b));
}

//S_uv
void calculateOverlapMatrix(BasisSet<STOFunction> a, arma::mat & Smatrix) {
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
    Smatrix.print("no rotation");
    int counter = 0;
    for (unsigned k = 0; k < a.myBasis.size(); k++) {
        for (unsigned l = 0; l < a.myBasis.size(); l++) {
            std::cout << k << l << std::endl;
            if (a.myBasis[k].id == a.myBasis[k - 1].id && a.myBasis[k].nlm.l == a.myBasis[k - 1].nlm.l && a.myBasis[k].id == a.myBasis[k + 1].id && a.myBasis[k].nlm.l == a.myBasis[k + 1].nlm.l && k > 1 && k < a.myBasis.size() - 1
                    && a.myBasis[l].id == a.myBasis[l - 1].id && a.myBasis[l].nlm.l == a.myBasis[l - 1].nlm.l && a.myBasis[l].id == a.myBasis[l + 1].id && a.myBasis[l].nlm.l == a.myBasis[l + 1].nlm.l && l > 1 && l < a.myBasis.size() - 1) {
                arma::mat rotationMatrix;
                findRotation(a.myBasis[k].x, a.myBasis[k].y, a.myBasis[k].z, a.myBasis[l].x, a.myBasis[l].y, a.myBasis[l].z, rotationMatrix);
                arma::mat sub(3, 3);
                sub = Smatrix.submat(k - 1, l - 1, k + 1, l + 1);
                sub = trans(rotationMatrix) * sub  * (rotationMatrix);
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        Smatrix(k - 1 + i, l - 1 + j) = sub(i, j);
                    }
                }
                counter++;
            }
            else if (a.myBasis[l].nlm.l == 0 && a.myBasis[k].id == a.myBasis[k - 1].id && a.myBasis[k].nlm.l == a.myBasis[k - 1].nlm.l && a.myBasis[k].id == a.myBasis[k + 1].id && a.myBasis[k].nlm.l == a.myBasis[k + 1].nlm.l && k > 1 && k < a.myBasis.size() - 1) {
                arma::mat rotationMatrix;
                findRotation(a.myBasis[k].x, a.myBasis[k].y, a.myBasis[k].z, a.myBasis[l].x, a.myBasis[l].y, a.myBasis[l].z, rotationMatrix);
                arma::mat sub(3, 1);
                sub = Smatrix.submat(k - 1, l, k + 1, l);
                sub = trans(rotationMatrix) * sub;
                for (int i = 0; i < 3; i++) {
                    Smatrix(k  + i - 1, l ) = sub(i, 0);
                }
                counter++;
            }
            else if (a.myBasis[k].nlm.l == 0 && a.myBasis[l].id == a.myBasis[l - 1].id && a.myBasis[l].nlm.l == a.myBasis[l - 1].nlm.l && a.myBasis[l].id == a.myBasis[l + 1].id && a.myBasis[l].nlm.l == a.myBasis[l + 1].nlm.l && l > 1 && l < a.myBasis.size() - 1) {
                arma::mat rotationMatrix;
                findRotation(a.myBasis[k].x, a.myBasis[k].y, a.myBasis[k].z, a.myBasis[l].x, a.myBasis[l].y, a.myBasis[l].z, rotationMatrix);
                arma::mat sub(1, 3);
                sub = (Smatrix.submat(k, l - 1, k, l + 1));
                sub = (sub * (rotationMatrix));
                for (int i = 0; i < 3; i++) {
                    Smatrix(k, l  + i - 1) = sub(0, i);
                }
                counter++;
            }
        }
    }
    //Smatrix(1, 5) = -Smatrix(1, 5);
    //ssSmatrix(5, 1) = -Smatrix(5, 1);
    std::cout << "counter:" << counter << std::endl;
}

} //namespace Semi
