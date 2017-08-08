#include "HuckelSolver.h"

namespace Semi {
void solveHuckelMatrixWithOverlap(arma::mat Hhuckel, arma::mat Shuckel, arma::mat Smatrix, Molecule m, std::vector<myOrbital> valenceOrbitalData, arma::mat &solMatrix) {
    arma::mat X;
    invSqrt(Shuckel, X);
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

    solMatrix = c_full;
}

void solveHuckelMatrix(arma::mat Hhuckel, arma::mat Smatrix, Molecule m, std::vector<myOrbital> valenceOrbitalData, arma::mat &solMatrix) {
    arma::mat eigvec;
    arma::vec eigval;
    eig_sym(eigval, eigvec, Hhuckel);

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
    solMatrix = eigvec;
}



} // namespace Semi
