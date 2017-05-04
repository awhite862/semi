#include "get_sto_3g.h"

namespace Semi {
void get_sto_3g(std::string &name, QNumber &qn, std::vector<double> &nVec, std::vector<double> &alphaVec) {
    if (name == "H" && qn.n == 1 && qn.l == 0) {
        alphaVec[0] = 3.4252509;   nVec[0] =  0.15432897;
        alphaVec[1] = 0.6239137;   nVec[1] =  0.53532814;
        alphaVec[2] = 0.1688554;   nVec[2] =  0.44463454;
    }
    else if (name == "He" &&  qn.n == 1 && qn.l == 0) {
        alphaVec[0] = 6.3624213;   nVec[0] =  0.15432897;
        alphaVec[1] = 1.1589230;   nVec[1] =  0.53532814;
        alphaVec[2] = 0.3136497;   nVec[2] =  0.44463454;
    }
    else if (name == "Li" &&  qn.n == 1 && qn.l == 0) {
        alphaVec[0] = 16.1195750;  nVec[0] =  0.15432897;
        alphaVec[1] =  2.9362007;  nVec[1] =  0.53532814;
        alphaVec[2] =  0.7946505;  nVec[2] =  0.44463454;
    }
    else if (name == "Li" &&  qn.n == 1 && qn.l == 0) {
        alphaVec[0] = 0.6362897;   nVec[0] = -0.09996723;
        alphaVec[1] = 0.1478601;   nVec[1] =  0.39951283;
        alphaVec[2] = 0.0480887;   nVec[2] =  0.70011547;
    }
    else if (name == "Li" &&  qn.n == 1 && qn.l == 1) {
        alphaVec[0] = 0.6362897;    nVec[0] =  0.15591627;
        alphaVec[1] = 0.1478601;    nVec[1] =  0.60768372;
        alphaVec[2] = 0.0480887;    nVec[2] =  0.39195739;
    }
    else if (name == "Be" && qn.n == 1 && qn.l == 0) {
        alphaVec[0] = 30.1678710;   nVec[0] =  0.15432897;
        alphaVec[1] =  5.4951153;   nVec[1] =  0.53532814;
        alphaVec[2] =  1.4871927;   nVec[2] =  0.44463454;
    }
    else if (name == "Be" && qn.n == 2 && qn.l == 1) {
        alphaVec[0] = 1.3148331;    nVec[0] =  0.15591627;
        alphaVec[1] = 0.3055389;    nVec[1] =  0.60768372;
        alphaVec[2] = 0.0993707;    nVec[2] =  0.39195739;
    }
    else if (name == "Be" && qn.n == 2 && qn.l == 0) {
        alphaVec[0] = 1.3148331;    nVec[0] = -0.09996723;
        alphaVec[1] = 0.3055389;    nVec[1] =  0.39951283;
        alphaVec[2] = 0.0993707;    nVec[2] =  0.70011547;
    }
    else if (name == "B" && qn.n == 1 && qn.l == 0) {
        alphaVec[0] = 48.7911130;   nVec[0] =  0.15432897;
        alphaVec[1] =  8.8873622;   nVec[1] =  0.53532814;
        alphaVec[2] =  2.4052670;   nVec[2] =  0.44463454;
    }
    else if (name == "B" && qn.n == 2 && qn.l == 1) {
        alphaVec[0] = 2.2369561;    nVec[0] = -0.09996723;
        alphaVec[1] = 0.5198205;    nVec[1] =  0.39951283;
        alphaVec[2] = 0.1690618;    nVec[2] =  0.70011547;
    }
    else if (name == "B" && qn.n == 2 && qn.l == 0) {
        alphaVec[0] = 2.2369561;    nVec[0] = -0.09996723;
        alphaVec[1] = 0.5198205;    nVec[1] =  0.39951283;
        alphaVec[2] = 0.1690618;    nVec[2] =  0.70011547;
    }
    else if (name == "C" && qn.n == 1 && qn.l == 0) {
        alphaVec[0] = 71.6168370;   nVec[0] =  0.15432897;
        alphaVec[1] = 13.0450960;   nVec[1] =  0.53532814;
        alphaVec[2] =  3.5305122;   nVec[2] =  0.44463454;
    }
    else if (name == "C" && qn.n == 2 && qn.l == 1) {
        alphaVec[0] = 2.9412494;    nVec[0] =  0.15591627;
        alphaVec[1] = 0.6834831;    nVec[1] =  0.60768372;
        alphaVec[2] = 0.2222899;    nVec[2] =  0.39195739;
    }
    else if (name == "C" && qn.n == 2 && qn.l == 0) {
        alphaVec[0] = 2.9412494;    nVec[0] = -0.09996723;
        alphaVec[1] = 0.6834831;    nVec[1] =  0.39951283;
        alphaVec[2] = 0.2222899;    nVec[2] =  0.70011547;
    }
    else if (name == "N" && qn.n == 1 && qn.l == 0) {
        alphaVec[0] = 99.1061690;   nVec[0] =  0.15432897;
        alphaVec[1] = 18.0523120;   nVec[1] =  0.53532814;
        alphaVec[2] =  4.8856602;   nVec[2] =  0.44463454;
    }
    else if (name == "N"  && qn.n == 2 && qn.l == 1) {
        alphaVec[0] = 3.7804559;    nVec[0] =  0.15591627;
        alphaVec[1] = 0.8784966;    nVec[1] =  0.60768372;
        alphaVec[2] = 0.2857144;    nVec[2] =  0.39195739;
    }
    else if (name == "N"  && qn.n == 2 && qn.l == 0) {
        alphaVec[0] = 3.7804559;    nVec[0] = -0.09996723;
        alphaVec[1] = 0.8784966;    nVec[1] =  0.39951283;
        alphaVec[2] = 0.2857144;    nVec[2] =  0.70011547;
    }
    else if (name == "O"  && qn.n == 1 && qn.l == 0) {
        alphaVec[0] = 130.7093200;  nVec[0] =  0.15432897;
        alphaVec[1] =  23.8088610;  nVec[1] =  0.53532814;
        alphaVec[2] =   6.4436083;  nVec[2] =  0.44463454;
    }
    else if (name == "O" && qn.n == 2 && qn.l == 1) {
        alphaVec[0] = 5.0331513;    nVec[0] =  0.15591627;
        alphaVec[1] = 1.1695961;    nVec[1] =  0.60768372;
        alphaVec[2] = 0.3803890;    nVec[2] =  0.39195739;
    }
    else if (name == "O" && qn.n == 2 && qn.l == 0) {
        alphaVec[0] = 5.0331513;    nVec[0] = -0.09996723;
        alphaVec[1] = 1.1695961;    nVec[1] =  0.39951283;
        alphaVec[2] = 0.3803890;    nVec[2] =  0.70011547;
    }
    else if (name == "F" && qn.n == 1 && qn.l == 0) {
        alphaVec[0] = 166.6791300;  nVec[0] =  0.15432897;
        alphaVec[1] =  30.3608120;  nVec[1] =  0.53532814;
        alphaVec[2] =   8.2168207;  nVec[2] =  0.44463454;
    }
    else if (name == "F" && qn.n == 2 && qn.l == 1) {
        alphaVec[0] = 6.4648032;    nVec[0] =  0.15591627;
        alphaVec[1] = 1.5022812;    nVec[1] =  0.60768372;
        alphaVec[2] = 0.4885885;    nVec[2] =  0.39195739;
    }
    else if (name == "F" && qn.n == 2 && qn.l == 0) {
        alphaVec[0] = 6.4648032;    nVec[0] = -0.09996723;
        alphaVec[1] = 1.5022812;    nVec[1] =  0.39951283;
        alphaVec[2] = 0.4885885;    nVec[2] =  0.70011547;
    }
    else if (name == "Ne" && qn.n == 1 && qn.l == 0) {
        alphaVec[0] = 207.0156100;  nVec[0] =  0.15432897;
        alphaVec[1] =  37.7081510;  nVec[1] =  0.53532814;
        alphaVec[2] =  10.2052970;  nVec[2] =  0.44463454;
    }
    else if (name == "Ne" && qn.n == 1 && qn.l == 0) {
        alphaVec[0] = 8.2463151;    nVec[0] = -0.09996723;
        alphaVec[1] = 1.9162662;    nVec[1] =  0.39951283;
        alphaVec[2] = 0.6232293;    nVec[2] =  0.70011547;
    }
    else {
        //throw std::runtime_error("Element out of range");
    }
}

} // namespace Semi