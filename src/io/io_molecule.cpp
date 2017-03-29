#include "io_molecule.h"
#include <iostream>
#include <sstream>
#include <iomanip>

namespace io {

molecular_input::molecular_input() { }

void molecular_input::print(std::ostream &os) const {

    for (size_t i = 0; i < vs.size(); i++) {
        os << std::setprecision(3);
        os << " " << vs[i] << " " << vx[i] << " " << vy[i] << " " << vz[i] << std::endl;
    }
}

bool molecular_input::is_empty() const {
    return vs.empty() && vx.empty() && vy.empty() && vz.empty();
}

void molecular_input::add(std::string &input) {
    std::stringstream ss;
    std::string s;
    double x,y,z;
    ss << input;
    ss >> s;
    ss >> x;
    ss >> y;
    ss >> z;
    vs.push_back(s);
    vx.push_back(x);
    vy.push_back(y);
    vz.push_back(z);
}

} // namespace io
