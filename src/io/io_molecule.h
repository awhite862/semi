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
    size_t size() const {return vs.size();}
    void add(std::string &input);
    std::string gets(size_t i) const  {return vs[i];}
    double getx(size_t i) const {return vx[i];}
    double gety(size_t i) const {return vy[i];}
    double getz(size_t i) const {return vz[i];}
};

} // namespace io

#endif // IO_MOLECULE_H
