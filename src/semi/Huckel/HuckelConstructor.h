#ifndef HUCKELCONSTRUCTOR_H
#define HUCKELCONSTRUCTOR_H
#include "HuckelMethod.h"
#include "semi/Structure/Atom.h"
#include "semi/Structure/Molecule.h"
#include "semi/semi_utils.h"
#include <armadillo>

namespace Semi {
/** \brief Class that solved Huckel Hamiltonain via diagonalization. **/

/** \brief Retrieves VOIE data from voie struct. **/
double fetchVOIE(double atom, std::string orbital, std::vector<voie> voieData);

/** \brief Retrieves all orbitals for a given element. **/
std::vector<std::string> fetchOrbitals(double atom);

/** \brief Retrieves all valence orbitals for a given element. **/
std::vector<std::string> fetchValenceOrbitals(double atom);

/** \brief Calculates k value from given parameters. **/
double kCalc(double k, double sigma, int i, int j, std::vector<myOrbital> valenceOrbitalData);

/** \brief Constructs a Huckel Hamiltonian from overlap matrix, xyz data, and parameters. **/
void constructHuckelHamiltonian(arma::mat Smatrix, double kValue, double sigmaValue, Molecule m,
                                arma::mat &Hhuckel, arma::mat &Shuckel, std::vector<myOrbital> &valenceOrbitalData);

} // namespace Semi

#endif // HUCKELCONSTRUCTOR_H
