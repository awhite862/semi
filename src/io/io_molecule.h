#ifndef IO_MOLECULE_H
#define IO_MOLECULE_H

#include <vector>
#include <string>

namespace io {

class molecular_input {
private:
    std::vector<std::string> vs;
    std::vector<double> vx, vy, vz;

public:
    molecular_input();
    void print(std::ostream &os) const;
    bool is_empty() const;
    void add(std::string &input);
};

} // namespace io

#endif // IO_MOLECULE_H
