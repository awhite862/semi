#define ARMA_DONT_USE_WRAPPER
#include <iostream>
#include <vector>
#include <armadillo>
#include <semi/semi_utils.h>
#include "huckel.h"
#include <tuple>
namespace Semi {


arma::mat invSqrt(arma::mat A) {
    arma::vec eigval;
    arma::mat eigvec;
    eig_sym(eigval, eigvec, inv(A));
    arma::mat eigvalmatrix = diagmat(eigval);
    eigvalmatrix = sqrt(eigvalmatrix);
    arma::mat ans = eigvec * eigvalmatrix * inv(eigvec);
    return ans;
}

double fetchVOIE(double atom, std::string orbital, std::vector<voie> voieData) {
    double H1;
    for (unsigned int k = 0; k < voieData.size(); k++) {
        if (voieData[k].atom == atom) {
            if (orbital.compare("1s") == 0) {
                H1 = voieData[k].oneS; break;
            }
            if (orbital.compare("2s") == 0) {
                H1 = voieData[k].twoS; break;
            }
            if (orbital.compare("2p") == 0) {
                H1 = voieData[k].twoP; break;
            }
            if (orbital.compare("3s") == 0) {
                H1 = voieData[k].threeS; break;
            }
            if (orbital.compare("3p") == 0) {
                H1 = voieData[k].threeP; break;
            }
        }
    }
    return H1;
}

std::vector<std::string> fetchOrbitals(double atom) {
    std::vector<std::string> v;
    if (atom <= 2) {
        const char* args[] = {"1s"};
        std::vector<std::string> v(args, args + 1);
        return v;
    }
    if (atom <= 10) {
        const char* args[] = {"1s", "2s", "2p", "2p", "2p"};
        std::vector<std::string> v(args, args + 5);
        return v;
    }
    if (atom <= 18) {
        const char* args[] = {"1s", "2s", "3s", "2p", "2p", "2p", "3p", "3p", "3p"};
        std::vector<std::string> v(args, args + 9);
        return v;
    }
    return v;
}

std::vector<std::string> fetchValenceOrbitals(double atom) {
    std::vector<std::string> v;
    if (atom <= 2) {
        const char* args[] = {"1s"};
        std::vector<std::string> v(args, args + 1);
        return v;
    }
    if (atom <= 10) {
        const char* args[] = {"2s", "2p", "2p", "2p"};
        std::vector<std::string> v(args, args + 4);
        return v;
    }
    if (atom <= 18) {
        const char* args[] = {"3s", "3p", "3p", "3p"};
        std::vector<std::string> v(args, args + 4);
        return v;
    }
    return v;
}

int numValence(int atom) {
    if (atom <= 2) {
        return atom;
    }
    if (atom <= 10) {
        return (atom - 2) % 8;;
    }
    if (atom <= 18) {
        return numValence(atom - 8);
    }
    return 0;
}

double kCalc(double k, double sigma, int i, int j, std::vector<myOrbital> valenceOrbitalData, std::vector<voie> voieData) {
    double num = (fetchVOIE(valenceOrbitalData[i].atom.charge, valenceOrbitalData[i].orbital, voieData) - fetchVOIE(valenceOrbitalData[j].atom.charge, valenceOrbitalData[j].orbital, voieData));
    double dem = (fetchVOIE(valenceOrbitalData[i].atom.charge, valenceOrbitalData[i].orbital, voieData) + fetchVOIE(valenceOrbitalData[j].atom.charge, valenceOrbitalData[j].orbital, voieData));
    double delta = num / dem;
    double dist = distance(valenceOrbitalData[i].atom.x, valenceOrbitalData[i].atom.y, valenceOrbitalData[i].atom.z, valenceOrbitalData[j].atom.x, valenceOrbitalData[j].atom.y, valenceOrbitalData[j].atom.z);
    double radii = zetaCalc(valenceOrbitalData[i].atom.charge) + zetaCalc(valenceOrbitalData[j].atom.charge);
    return (1.0 + k + pow(delta, 2) - pow(delta, 4)) * exp(sigma * (dist - radii));
}

void get_VOIE(unsigned data, double &v1, double &v2, double &v3, double &v4, double &v5) {

    switch (data) {
    case 1:   v1 = 13.6; v2 = 0;    v3 = 0;    v4 = 0;    v5 = 0;       break;
    case 2:   v1 = 24.5; v2 = 0;    v3 = 0;    v4 = 0;    v5 = 0;       break;
    case 3:   v1 = 0;    v2 = 5.45; v3 = 3.5;  v4 = 0;    v5 = 0;       break;
    case 4:   v1 = 0;    v2 = 9.30; v3 = 6.0;  v4 = 0;    v5 = 0;       break;
    case 5:   v1 = 0;    v2 = 14.0; v3 = 8.30; v4 = 0;    v5 = 0;       break;
    case 6:   v1 = 0;    v2 = 19.5; v3 = 10.7; v4 = 0;    v5 = 0;       break;
    case 7:   v1 = 0;    v2 = 25.5; v3 = 13.1; v4 = 0;    v5 = 0;       break;
    case 8:   v1 = 0;    v2 = 32.3; v3 = 15.9; v4 = 0;    v5 = 0;       break;
    case 9:   v1 = 0;    v2 = 40.4; v3 = 18.7; v4 = 0;    v5 = 0;       break;
    case 10:  v1 = 0;    v2 = 48.5; v3 = 21.5; v4 = 0;    v5 = 0;       break;
    case 11:  v1 = 0;    v2 = 0;    v3 = 0;    v4 = 5.21; v5 = 0;       break;
    case 12:  v1 = 0;    v2 = 0;    v3 = 0;    v4 = 7.68; v5 = 0;       break;
    case 13:  v1 = 0;    v2 = 0;    v3 = 0;    v4 = 11.3; v5 = 5.95;    break;
    case 14:  v1 = 0;    v2 = 0;    v3 = 0;    v4 = 15.0; v5 = 7.81;    break;
    case 15:  v1 = 0;    v2 = 0;    v3 = 0;    v4 = 18.7; v5 = 10.2;    break;
    case 16:  v1 = 0;    v2 = 0;    v3 = 0;    v4 = 20.7; v5 = 11.7;    break;
    case 17:  v1 = 0;    v2 = 0;    v3 = 0;    v4 = 25.3; v5 = 13.8;    break;
    }
}

#define EV_TO_HARTREE 27.21138602

std::tuple <arma::mat, arma::mat, std::vector<myOrbital>> constructHuckelHamiltonian(arma::mat Smatrix, double kValue, double sigmaValue, Molecule m) {
    int i = 0;
    std::vector<voie> voieData;
    std::vector<myOrbital> allOrbitalData;
    std::vector<myOrbital> valenceOrbitalData;
    for (size_t i = 0; i < 18; i++) {
        double scale = EV_TO_HARTREE;
        double val1, val2, val3, val4, val5;
        get_VOIE(i, val1, val2, val3, val4, val5);
        voieData.push_back(voie());
        voieData[i].atom = i; voieData[i].oneS = val1 / scale; voieData[i].twoS = val2 / scale;
        voieData[i].twoP = val3 / scale; voieData[i].threeS = val4 / scale; voieData[i].threeP = val5 / scale;
    }

    i = 0;
    for (unsigned int k = 0; k < m.myMolecule.size(); k++) {
        std::vector<std::string> orbital_no_atom = fetchOrbitals(m.myMolecule[k].charge);
        for (unsigned int j = 0; j < orbital_no_atom.size(); j++) {
            allOrbitalData.push_back(myOrbital());
            allOrbitalData[i].atom = Atom(m.myMolecule[k].x, m.myMolecule[k].y, m.myMolecule[k].z, m.myMolecule[k].charge, k);
            allOrbitalData[i].orbital = orbital_no_atom[j];
            i++;
        }
    }

    i = 0;
    for (unsigned int k = 0; k < m.myMolecule.size(); k++) {
        std::vector<std::string> valence_orbital_no_atom = fetchValenceOrbitals(m.myMolecule[k].charge);
        for (unsigned int j = 0; j < valence_orbital_no_atom.size(); j++) {
            valenceOrbitalData.push_back(myOrbital());
            valenceOrbitalData[i].atom = Atom(m.myMolecule[k].x, m.myMolecule[k].y, m.myMolecule[k].z, m.myMolecule[k].charge, k);
            valenceOrbitalData[i].orbital = valence_orbital_no_atom[j];
            i++;
        }
    }

    std::vector<double> keep;
    unsigned int size = sqrt(Smatrix.size());
    if (valenceOrbitalData.size() != size) {
        for (unsigned int i = 0; i < size; i++) {
            std::vector<std::string> v = fetchValenceOrbitals(allOrbitalData[i].atom.charge);
            if ((find(v.begin(), v.end(), allOrbitalData[i].orbital) != v.end())) {
                keep.push_back(i);
            }
        }
    }
    else {
        for (unsigned int i = 0; i < size; i++) {
            keep.push_back(i);
        }
    }

    arma::Col<arma::uword> vec(keep.size());
    for (unsigned int k = 0; k < keep.size(); k++) {
        vec[k] = keep[k];
    }

    arma::mat Shuckel(vec.size(), vec.size());
    Shuckel = Smatrix.submat(vec, vec);

    arma::mat Hhuckel(valenceOrbitalData.size(), valenceOrbitalData.size());
    for (unsigned int i = 0; i < keep.size(); i++) {
        for (unsigned int j = 0; j < keep.size(); j++) {
            if (i != j) {
                Hhuckel(i, j) = -0.5 * kCalc(kValue, sigmaValue, i, j, valenceOrbitalData, voieData) * Shuckel(i, j) * (fetchVOIE(valenceOrbitalData[i].atom.charge, valenceOrbitalData[i].orbital, voieData) + fetchVOIE(valenceOrbitalData[j].atom.charge, valenceOrbitalData[j].orbital, voieData));
            }
            else {
                Hhuckel(i, i) = -fetchVOIE(valenceOrbitalData[i].atom.charge, valenceOrbitalData[i].orbital, voieData);
            }
        }
    }
    Hhuckel.print("Here is H");

    return std::make_tuple(Hhuckel, Shuckel, valenceOrbitalData);
}

arma::mat solveHuckelMatrix(arma::mat Hhuckel, arma::mat Shuckel, arma::mat Smatrix, Molecule m, std::vector<myOrbital> valenceOrbitalData) {
    arma::mat X = invSqrt(Shuckel);
    arma::mat Hprime = X * Hhuckel * X;
    arma::mat eigvec;
    arma::vec eigval;
    eig_sym(eigval, eigvec, Hprime);
    eigvec = X * eigvec;

    unsigned int numpi = 0;
    for (unsigned int i = 0; i < m.myMolecule.size(); i++) {
        numpi += numValence(m.myMolecule[i].charge);
    }

    unsigned int num_orbitals = numpi / 2;

    arma::Col<arma::uword> retain(num_orbitals);
    for (unsigned int i = 0; i < num_orbitals; i++) {
        retain[i] = i;
    }

    arma::mat c_v(valenceOrbitalData.size(), num_orbitals);
    c_v = eigvec.cols(0, num_orbitals - 1);

    unsigned int size = sqrt(Smatrix.size());
    arma::mat c_full(size, size - valenceOrbitalData.size() + num_orbitals);
    for (unsigned int k = 0; k < size; k++) {
        for (unsigned int i = 0; i < size - valenceOrbitalData.size() + num_orbitals; i++) {
            c_full(k, i) = 0;
        }
    }
    for (unsigned int k = 0; k < size - valenceOrbitalData.size(); k++) {
        c_full(k, k) = 1;
    }
    for (unsigned int k = 0; k < valenceOrbitalData.size(); k++) {
        for (unsigned int i = 0; i < num_orbitals; i++) {
            c_full(k + size - valenceOrbitalData.size(), i + size - valenceOrbitalData.size()) = c_v(k, i);
        }
    }


    std::cout << "Orbital Energies" << std::endl;
    for (int k = 0; k < eigval.size(); k++) {
        if (eigval(k) > 0)
            std::cout << "Energy = " << "+" << eigval(k) << std::endl;
        else
            std::cout << "Energy = " << "-" << (-1.00 * eigval(k)) << std::endl;
    }

    double energy = 0;
    for (int k = 0; k < num_orbitals; k++) {
        energy += 2.0 * eigval(k);
    }
    std::cout << "Total Energy: " << energy << std::endl;
}

arma::mat constructCBasis(arma::mat c_full, arma::mat Smatrix) {
    arma::mat y = trans(c_full) * Smatrix * c_full;
    arma::mat c_basis = c_full * (invSqrt(y));
    arma::mat id = trans(c_basis) * Smatrix * c_basis;
    return c_basis;
}
arma::mat calculateHuckel(arma::mat Smatrix, double kValue, double sigmaValue, Molecule m, std::string args) {
    auto t = constructHuckelHamiltonian(Smatrix, kValue, sigmaValue, m);
    arma::mat a= std::get<0>(t);
    solveHuckelMatrix(std::get<0>(t), std::get<1>(t), Smatrix, m, std::get<2>(t));
}

}