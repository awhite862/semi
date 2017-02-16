#ifndef BASIS_H
#define BASIS_H

#include <string>
#include <sstream>

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
	Basis(double zetaValue, int nValue, int lValue, int mValue, double xValue, double yValue, double zValue);

	//get string represeting orbital type
	std::string getOrbitalType();



};
}

#endif //BASIS_H
