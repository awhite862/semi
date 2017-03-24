/** \brief Class that represents approxiamtes initial guess via Huckel theory. **/
#ifndef IntegralEvaluator_H
#define IntegralEvaluator_H
#include <armadillo>
#include <vector>

namespace Semi {

/** \brief Struct that stores voie data. **/
struct voie;

/** \brief Struct that stores xyz and element info. **/
struct xyz;

/** \brief Calculates the inverse square root of a matrix. **/
struct myOrbital;

/** \brief Calculates the inverse square root of a matrix. **/
arma::mat invSqrt(arma::mat A);

/** \brief Retrieves VOIE data from voie struct. **/
double fetchVOIE(std::string atom1, std::string orbital1);

/** \brief Retrieves all orbitals for a given element. **/
vector<std::string> fetchOrbitals(double atom);

/** \brief Retrieves all valence orbitals for a given element. **/
vector<string> fetchValenceOrbitals(double atom);

/** \brief Calculates number of valence electrons for given atom. **/
int numValence(int num);

/** \brief Calculates k value from given parameters. **/
double kCalc(double k, double sigma, int i, int j);

/** \brief Huckel approximation of initial guess given overlap matrix, xyz data and parameters. **/
mat huckel(mat Smatrix, double kValue, double sigmaValue, vector<xyz> xyzData);

/** \brief Main. **/
int main();

} //namespace Semi
#endif //huckel.h