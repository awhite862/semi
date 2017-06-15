set(SRC
    Basis/GTOFunction.cpp
    Basis/CGTOFunction.cpp
    Basis/STOFunction.cpp
    Basis/get_sto_3g.cpp
    Huckel/HuckelConstructor.cpp
    Huckel/HuckelSolver.cpp
    Huckel/HuckelMethod.cpp
    Huckel/HuckelBasis.cpp
    Integral/Cndo.cpp
    Integral/IntegralEvaluator.cpp
    Structure/Atom.cpp
    Structure/QNumber.cpp
    Structure/Molecule.cpp
    run_huckel.cpp
    semi_utils.cpp
    semi_api.cpp
)

add_library(libsemi ${SRC})
include_directories(${SEMI_DIR}/src)
