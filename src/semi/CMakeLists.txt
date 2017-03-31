set(SRC
    Atom.cpp
    Basis.cpp
    BasisSet.cpp
    Molecule.cpp
    Integral/IntegralEvaluator.cpp
    Integral/Cndo.cpp
    Huckel/huckel.cpp
    semi_api.cpp
)

add_library(libsemi ${SRC})
include_directories(${SEMI_DIR}/src)
