#include "semi_utils.h"
#include "semi/Atom.h"

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
    return "null";
}

int factorial(int n) {
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

double delta(double i, double j) {
    return i == j ? 1 : 0;
}

} // namespace Semi