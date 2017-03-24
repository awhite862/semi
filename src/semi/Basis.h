/** \brief Class that represents a basis. **/
#ifndef BASIS_H
#define BASIS_H

namespace Semi {
class Basis {
public:
/** \brief Zeta value of basis. **/
	double zeta;

/** \brief n quantum number of basis. **/
	int n;

/** \brief l quantum number of basis. **/
	int l;

/** \brief m quantum number of basis. **/
	int m;

/** \brief x value of basis. **/
	double x;

/** \brief y value of basis. **/
	double y;

/** \brief z value of basis. **/
	double z;

/** \brief Constructor. **/
	Basis(double _zeta, int _n, int _l, int _m, double _x, double _y, double _z);
};

} // namespace Semi

#endif //BASIS_H
