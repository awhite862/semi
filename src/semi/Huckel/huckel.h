/** \brief Class that represents approxiamtes initial guess via Huckel theory. **/
#ifndef HUCKEL_H
#define HUCKEL_H
#include <armadillo>
#include <vector>
#include "semi/Molecule.h"
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
    Atom atom;
    std::string orbital;
};

/** \brief Calculates the inverse square root of a matrix. **/
arma::mat invSqrt(arma::mat A);

/** \brief Retrieves VOIE data from voie struct. **/
double fetchVOIE(double atom, std::string orbital, std::vector<voie> voieData);

/** \brief Retrieves all orbitals for a given element. **/
std::vector<std::string> fetchOrbitals(double atom);

/** \brief Retrieves all valence orbitals for a given element. **/
std::vector<std::string> fetchValenceOrbitals(double atom);

/** \brief Calculates number of valence electrons for given atom. **/
int numValence(int num);

/** \brief Calculates k value from given parameters. **/
double kCalc(double k, double sigma, int i, int j, std::vector<myOrbital> valenceOrbitalData);

/** \brief Huckel approximation of initial guess given overlap matrix, xyz data and parameters.
 **  args controls the output matrix.
**/
arma::mat calculateHuckel(arma::mat Smatrix, double kValue, double sigmaValue, Molecule m, std::string args);

} //namespace Semi
#endif //huckel.h
