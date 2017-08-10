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
    run_cndo.cpp
    semi_utils.cpp
    semi_api.cpp
    cndo2.cpp
    scf.cpp
    run_cndo.cpp
    run_scf.cpp 
)

add_library(libsemi ${SRC})
include_directories(${SEMI_DIR}/src)
