#ifndef QNUMBER_H
#define QNUMBER_H

namespace Semi {
/** \brief Class that represents an set of quantum numbers. **/
class QNumber {
public:
    int n; ///!< n quantum number
    int l; ///!< l quantum number
    int m; ///!< m quantum number

public:
    /** \brief Constructor.
        \param _n n.
        \param _l l.
        \param _m m.
     **/
    QNumber(int _n, int _l, int _m);

    /** \brief Default constructor.
     **/
    QNumber();
};

} // namespace Semi

#endif // QNUMBER_H
