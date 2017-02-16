set(SRC
    Atom.cpp
    Basis.cpp
    BasisSet.cpp
    Molecule.cpp
)

add_library(libsemi ${SRC})
include_directories(${SEMI_DIR}/src)
