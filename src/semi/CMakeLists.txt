set(SRC
    Atom.cpp
    STOBasis.cpp
    GTOBasis.cpp
    CGTOBasis.cpp
    QNumber.cpp
    get_sto_3g.cpp
    Molecule.cpp
    Integral/IntegralEvaluator.cpp
    Integral/Cndo.cpp
    Huckel/huckel.cpp
    run_huckel.cpp
    semi_utils.cpp
    semi_api.cpp
)

add_library(libsemi ${SRC})
include_directories(${SEMI_DIR}/src)
