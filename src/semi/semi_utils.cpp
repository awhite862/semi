#include "semi_utils.h"
#include <cstdlib>
#include <stdexcept>
namespace Semi {
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

double zetaCalc(double charge){
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

int factorial(int n) {
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}
 
int doubleFactorial(int n){
    return (n <= 1) ? 1 : factorial(n - 2) * n;
}

double distance (double x1, double y1, double z1, double x2, double y2, double z2) {
    return sqrt((pow((x1 - x2), 2) + pow((y1 - y2), 2) + pow((z1 - z2), 2)));
}

double delta(double i, double j) {
    return i == j ? 1 : 0;
}

} // namespace Semi