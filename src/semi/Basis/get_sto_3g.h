#ifndef GET_STO_3G_H
#define GET_STO_3G_H

#include <string>
#include <vector>
#include "semi/Structure/QNumber.h"

namespace Semi {

/** \brief Fetches STO-3g constraction coefficients and alpha terms for given element and quantum number.
        \param 
        \param _a x power.
        \param _b y power.
        \param _c z power.
        \param _r col vector representing vector.
        \param _elem atomic number.
        \param _nlm Quantum numbers nlm.
     **/

void get_sto_3g(std::string &name, QNumber &qn, std::vector<double> &nVec, std::vector<double> &alphaVec);

} // namespace Semi

#endif // GET_STO_3G_H