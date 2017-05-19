#include <string>
#include <vector>
#include <semi/QNumber.h>

namespace Semi {

/** \brief Fetches STO-3g constraction coefficients and alpha terms for given element and quantum number. **/
void get_sto_3g(std::string &name, QNumber &qn, std::vector<double> &nVec, std::vector<double> &alphaVec);
}
