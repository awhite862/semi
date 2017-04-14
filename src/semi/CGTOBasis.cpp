#include "CGTOBasis.h"
#include <armadillo>
#include <iostream>
#include <map>
#include <math.h>
#include "GTOBasis.h"
#include "QNumber.h"

namespace Semi {
CGTOBasis::CGTOBasis(double _a, double _b, double _c, arma::colvec _r, double _elem, QNumber _nlm) {
    a = _a;
    b = _b;
    c = _c;
    r = _r;
    elem = _elem;
    nlm = _nlm;

    l = a + b + c;

    std::ifstream inputFile("sto-3g.txt");
    std::string line;
    int counter = 0;
    std::string key;
    std::map<std::string, std::vector<coeffs> > stoMap;
    std::vector<coeffs> temp;
    while (std::getline(inputFile, line)) {
        //std::cout << line << std::endl;
        if (counter == 0) {
            std::stringstream ss(line);
            std::string var1, var2;
            ss >> var1 >> var2;
            key = var1 + var2;
            ///std::cout << key << std::endl;
        }
        else if (counter % 4 == 0) {
            stoMap[key] = temp;
            temp.clear();
            std::stringstream ss(line);
            std::string var1, var2;
            ss >> var1 >> var2;
            key = var1 + var2;
            ///std::cout << key << std::endl;
        }
        else {
            std::stringstream ss(line);
            double var1, var2, var3;
            ss >> var1 >> var2 >> var3;
            ///std::cout << var1 << " " << var2 << " " << var3 << std::endl;
            coeffs s = coeffs();
            s.a = var1;
            s.cs = var2;
            s.cp = var3;
            temp.push_back(s);
        }
        counter++;
    }
    stoMap[key] = temp;

    //std::map<std::string, std::vector<coeffs> > stomap = preload();
    int z = (int) (elem + 0.1);
    std::string name;
    switch ((int) z) {
    case 1: name = "H"; break;
    case 2: name = "He"; break;
    case 3: name = "Li"; break;
    case 4: name = "Be"; break;
    case 5: name = "B"; break;
    case 6: name = "C"; break;
    case 7: name = "N"; break;
    case 8: name = "O"; break;
    case 9: name = "F"; break;
    case 10: name = "Ne"; break;
    }

    std::ostringstream oss;
    oss << nlm.n;
    std::string keys = name + oss.str() + "S";
    if (nlm.l != 0) {
        keys += "P";
    }
    double n;
    std::cout << keys << std::endl;
    std::vector<coeffs> coeffVector = stoMap[keys];
    for (int k = 0; k < 3; k++) {
        alphaVec.push_back(coeffVector[k].a);
        n=1;
        //n =  M_PI / (2 * coeffVector[k].a) * pow((factorial(factorial(2.0 * a - 1)) * factorial(factorial(2.0 * b - 1)) * factorial(factorial(2.0 * c - 1))) / (pow(2, 2 * l) * pow(coeffVector[k].a, l)), -1.0 / 2.0);
        if (nlm.l = 0) {
            nVec.push_back(coeffVector[k].cs * n);
        }
        else {
            nVec.push_back(coeffVector[k].cp * n);
        }
    }
    //std::cout << coeffVector[0].a << " " <<  coeffVector[0].cs  << " " <<  coeffVector[0].cs << std::endl;
}



} //namespace Semi
