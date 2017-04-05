set(SRC
    Atom.cpp
    STOBasis.cpp
    GTOBasis.cpp
    BasisSet.cpp
    Molecule.cpp
    Integral/IntegralEvaluator.cpp
    Integral/Cndo.cpp
    Huckel/huckel.cpp
)

add_library(libsemi ${SRC})
include_directories(${SEMI_DIR}/src)
