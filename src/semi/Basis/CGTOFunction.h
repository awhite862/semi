#ifndef CGTOFUNCTION_H
#define CGTOFUNCTION_H

#include <semi/QNumber.h>
#include <armadillo>

namespace Semi {
/** \brief Struct used to read in sto-3g data from EMSL . **/
struct coeffs {
    double a;
    double cs;
    double cp;
};

/** \brief Class that represents an CGTO Basis Function. **/
class CGTOFunction {
public:
    double a; ///!< x power
    double b; ///!< y power
    double c; ///!< z power
    std::vector<double> nVec; ///!< vector for normalization
    std::vector<double> alphaVec; ///!< vector for alpha values
    arma::colvec r; ///!< col vector representing location
    QNumber nlm; ///!< Quantum numbers nlm
    double elem; ///!< atomic number

public:
    /** \brief Constructor with all values passed in explicitly.
        \param _alphavec vector of alpha values
        \param _a x power.
        \param _b y power.
        \param _c z power.
        \param _r col vector representing vector.
        \param _elem atomic number.
        \param _nlm Quantum numbers nlm.
     **/
    CGTOFunction(QNumber _nlm, double _a, double _b, double _c, arma::colvec _r, double _elem, std::vector<double> &_nVec, std::vector<double> &_alphaVec) :
        nlm(_nlm), a(_a), b(_b), c(_c), r(_r), elem(_elem), nVec(_nVec), alphaVec(_alphaVec) { }

    /** \brief Constructor that requires minimum number of values to initiate cgto basis and calculates the rest.
        \param _a x power.
        \param _b y power.
        \param _c z power.
        \param _r col vector representing vector.
        \param _elem atomic number.
        \param _nlm Quantum numbers nlm.
     **/
    CGTOBasis(QNumber _nlm, double _a, double _b, double _c, arma::colvec _r, double _elem);
};

} // namespace Semi

#endif // GTOBASIS_H
