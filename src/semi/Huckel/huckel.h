/** \brief Class that represents approxiamtes initial guess via Huckel theory. **/
#ifndef HUCKEL_H
#define HUCKEL_H
#include <armadillo>
#include <vector>

namespace Semi {

/** \brief Struct that stores xyz and element info. **/
struct xyz {
	double atom;
	double x;
	double y;
	double z;
};

/** \brief Calculates the inverse square root of a matrix. **/
arma::mat invSqrt(arma::mat A);

/** \brief Retrieves VOIE data from voie struct. **/
double fetchVOIE(std::string atom1, std::string orbital1);

/** \brief Retrieves all orbitals for a given element. **/
std::vector<std::string> fetchOrbitals(double atom);

/** \brief Retrieves all valence orbitals for a given element. **/
std::vector<std::string> fetchValenceOrbitals(double atom);

/** \brief Calculates number of valence electrons for given atom. **/
int numValence(int num);

/** \brief Calculates k value from given parameters. **/
double kCalc(double k, double sigma, int i, int j);

/** \brief Huckel approximation of initial guess given overlap matrix, xyz data and parameters.
 **  args controls the output matrix.
**/
arma::mat calculateHuckel(arma::mat Smatrix, double kValue, double sigmaValue, std::vector<xyz> xyzData, std::string args);

} //namespace Semi
#endif //huckel.h
