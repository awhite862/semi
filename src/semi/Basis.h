#ifndef BASIS_H
#define BASIS_H

#include <string>

namespace Semi {
class Basis {
public:
	double zeta;
	int n;
	int l;
	int m;
	double x;
	double y;
	double z;

	//constructor
	Basis(double _zeta, int _n, int _l, int _m, double _x, double _y, double _z);

	//get string represeting orbital type
	std::string getOrbitalType();
};

} // namespace Semi

#endif //BASIS_H
