#define ARMA_DONT_USE_WRAPPER
#include <iostream>
#include <vector>
#include <armadillo>
#include "huckel.h"
namespace Semi {

struct voie {
    double atom;
    double oneS;
    double twoS;
    double twoP;
    double threeS;
    double threeP;
};

struct myOrbital {
    double atom;
    std::string orbital;
};

std::vector<double> atoms;
std::vector<voie> voieData;
std::vector<myOrbital> allOrbitalData;
std::vector<myOrbital> valenceOrbitalData;

arma::mat invSqrt(arma::mat A) {
    arma::vec eigval;
    arma::mat eigvec;
    eig_sym(eigval, eigvec, inv(A));
    arma::mat eigvalmatrix = diagmat(eigval);
    eigvalmatrix = sqrt(eigvalmatrix);
    arma::mat ans = eigvec * eigvalmatrix * inv(eigvec);
    return ans;
}

double fetchVOIE(double atom, std::string orbital) {
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

double kCalc(double k, double sigma, int i, int j) {
    double num = (fetchVOIE(valenceOrbitalData[i].atom, valenceOrbitalData[i].orbital) - fetchVOIE(valenceOrbitalData[j].atom, valenceOrbitalData[j].orbital));
    double dem = (fetchVOIE(valenceOrbitalData[i].atom, valenceOrbitalData[i].orbital) + fetchVOIE(valenceOrbitalData[j].atom, valenceOrbitalData[j].orbital));
    double delta = num / dem;
    return 1.0 + (k + pow(delta, 2) - pow(delta, 4) * k);
}

arma::mat calculateHuckel(arma::mat Smatrix, double kValue, double sigmaValue, std::vector<xyz> xyzData, std::string args) {
    std::ifstream fin;
    int i = 0;
    std::string label;
    std::string line;

    fin.open("VOIE.txt");
    while (std::getline(fin, line)) {
        std::stringstream linestream(line);
        double data, val1, val2, val3, val4, val5;
        linestream >> data >> val1 >> val2 >> val3 >> val4 >> val5;
        voieData.push_back(voie());
        voieData[i].atom = data; voieData[i].oneS = val1; voieData[i].twoS = val2;
        voieData[i].twoP = val3; voieData[i].threeS = val4; voieData[i].threeP = val5;
        i++;
    }
    fin.close();

    for (unsigned int k = 0; k < xyzData.size(); k++) {
        atoms.push_back(xyzData[k].atom);
    }

    i = 0;
    for (unsigned int k = 0; k < atoms.size(); k++) {
        std::vector<std::string> orbital_no_atom = fetchOrbitals(atoms[k]);
        for (unsigned int j = 0; j < orbital_no_atom.size(); j++) {
            allOrbitalData.push_back(myOrbital());
            allOrbitalData[i].atom = atoms[k];
            allOrbitalData[i].orbital = orbital_no_atom[j];
            i++;
        }
    }

    i = 0;
    for (unsigned int k = 0; k < atoms.size(); k++) {
        std::vector<std::string> valence_orbital_no_atom = fetchValenceOrbitals(atoms[k]);
        for (unsigned int j = 0; j < valence_orbital_no_atom.size(); j++) {
            valenceOrbitalData.push_back(myOrbital());
            valenceOrbitalData[i].atom = atoms[k];
            valenceOrbitalData[i].orbital = valence_orbital_no_atom[j];
            i++;
        }
    }

    std::vector<double> keep;
    unsigned int size = sqrt(Smatrix.size());
    for (unsigned int i = 0; i < size; i++) {
        std::vector<std::string> v = fetchValenceOrbitals(allOrbitalData[i].atom);
        if ((find(v.begin(), v.end(), allOrbitalData[i].orbital) != v.end())) {
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
                Hhuckel(i, j) = -0.5 * kCalc(kValue, sigmaValue, i, j) * Shuckel(i, j) * (fetchVOIE(valenceOrbitalData[i].atom, valenceOrbitalData[i].orbital) + fetchVOIE(valenceOrbitalData[j].atom, valenceOrbitalData[j].orbital));
            }
            else {
                Hhuckel(i, i) = -fetchVOIE(valenceOrbitalData[i].atom, valenceOrbitalData[i].orbital);
            }
        }
    }

    arma::mat eigvec;
    arma::vec eigval;
    eig_sym(eigval, eigvec, Hhuckel);

    unsigned int numpi = 0;
    for (unsigned int i = 0; i < atoms.size(); i++) {
        numpi += numValence(atoms[i]);
    }

    unsigned int num_orbitals = numpi / 2;

    arma::Col<arma::uword> retain(num_orbitals);
    for (unsigned int i = 0; i < num_orbitals; i++) {
        retain[i] = i;
    }

    arma::mat c_v(valenceOrbitalData.size(), num_orbitals);
    c_v = eigvec.cols(0, num_orbitals - 1);

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

    arma::mat y = trans(c_full) * Smatrix * c_full;
    arma::mat c_basis = c_full * (invSqrt(y));
    arma::mat id = trans(c_basis) * Smatrix * c_basis;

    /*  std::cout << "eigenvalue" << std::endl;
        eigval.print();
        std::cout << "eigenvalue size " << eigval.n_elem << std::endl;
        std::cout << "eigenstd::vector" << std::endl;
        eigvec.print();
        std::cout << "eigenstd::vector size " << eigvec.n_elem << std::endl;
        std::cout << "numpi" << numpi << std::endl;
        std::cout << "number of valence" << std::endl;
        std::cout << valenceOrbitalData.size() << std::endl;
        std::cout << "number of filled orbitals" << std::endl;
        std::cout << num_orbitals << std::endl;
        std::cout << "c_v size " << c_v.n_elem << std::endl;
        c_v.print();
        std::cout << "c_full" << c_full.n_elem << std::endl;
        c_full.print();
        std::cout << "y" << std::endl;
        y.print();
        std::cout << "basis matrix" << std::endl;
        c_basis.print();
        std::cout << "identity matrix check" << std::endl;
        id.print();
    */

    if (!args.compare("c_v")) {
        return c_v;
    }
    else if (!args.compare("c_full")) {
        return c_full;
    }
    else if (!args.compare("c_basis")) {
        return c_basis;
    }
} l
}