#include "semi_utils.h"
#include <cstdlib>
#include <stdexcept>
namespace Semi {

//basic math functions
int factorial(int n) {
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

int doubleFactorial(int n) {
    return (n <= 1) ? 1 : factorial(n - 2) * n;
}

double delta(double i, double j) {
    return i == j ? 1 : 0;
}

double distance (double x1, double y1, double z1, double x2, double y2, double z2) {
    return sqrt((pow((x1 - x2), 2) + pow((y1 - y2), 2) + pow((z1 - z2), 2)));
}

//armadillo extra functionality
void invSqrt(arma::mat A, arma::mat &sol) {
    arma::vec eigval;
    arma::mat eigvec;
    eig_sym(eigval, eigvec, inv(A));
    arma::mat eigvalmatrix = diagmat(eigval);
    eigvalmatrix = sqrt(eigvalmatrix);
    arma::mat ans = eigvec * eigvalmatrix * inv(eigvec);
    sol = ans;
}

//misc chemistry related functions
std::string getElement(double charge) {
    int temp = (int) (charge + 0.1);
    switch ((int) temp) {
    case 1: return "H";
    case 2: return "He";
    case 3: return "Li";
    case 4: return "Be";
    case 5: return "B";
    case 6: return "C";
    case 7: return "N";
    case 8: return "O";
    case 9: return "F";
    case 10: return "Ne";
    }
    throw std::runtime_error("Element out of range");
}

double zetaCalc(double charge) {
    int temp = (int) (charge + 0.1);
    switch ((int) temp) {
    case 1: return 1.0;
    case 2: return 1.7;
    case 3: return 0.65;
    case 4: return 0.975;
    case 5: return 1.3;
    case 6: return 1.625;
    case 7: return 1.95;
    case 8: return 2.275;
    case 9: return 2.6;
    case 10: return 2.925;
    }
    throw std::runtime_error("Element out of range: " + std::to_string(temp));
}

int numValence(int atom) {
    if (atom <= 2) {
        return atom;
    }
    else if (atom <= 10) {
        return (atom - 2);
    }
    else {
        return 100;
    }
}

double getIE(double a) {
    double h = 27.21138602;
    switch ((int)a) {
    case 1: return 13.5984 / h;
    case 2: return 24.5874 / h;
    case 3: return 5.3917 / h;
    case 4: return 9.3227 / h;
    case 5: return 8.298 / h;
    case 6: return 11.2603 / h;
    case 7: return 14.5341 / h;
    case 8: return 13.6181 / h;
    case 9: return 17.4228 / h;
    case 10: return 0;
    }
    return 0;
}

double getEA(double a) {
    double h = 27.21138602;
    switch ((int)a) {
    case 1: return 0.754 / h;
    case 2: return -19.7 / h;
    case 3: return 0.618 / h;
    case 4: return -2.4 / h;
    case 5: return 0.279 / h;
    case 6: return 1.262 / h;
    case 7: return -1.4 / h;
    case 8: return 1.461 / h;
    case 9: return 3.401 / h;
    case 10: return 0;
    }
    return 0;
}

// int numPi(BasisSet<STOFunction> a) {
//     unsigned int sum = 0;
//     std::vector<int> ids;
//     for (int k = 0; k < a.myBasis.size(); k++) {
//         if (std::find(ids.begin(), ids.end(), a.myBasis[k].id) == ids.end()) {
//             ids.push_back(a.myBasis[k].id);
//             numpi += numValence(a.myBasis[k].charge);
//         }
//     }
//     return sum;
// }



} // namespace Semi