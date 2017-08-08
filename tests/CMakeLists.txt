set(TESTS 
    arma_test
    io_utils_test
    framework_test
    semi_test
    sto_test
    cgto_test
    sto_vs_cgto_test
    huckel_test
    cndo_test
)

#find_package( BLAS REQUIRED )
#find_package( LAPACK REQUIRED )

include_directories(${SEMI_DIR}/src)
include_directories(${LAPACK_INCLUDE_DIR}${BLAS_INCLUDE_DIR}${ARMADILLO_INCLUDE_DIRS})

foreach(TEST ${TESTS})
    add_executable(${TEST} ${TEST}.cpp)   
    target_link_libraries(${TEST})
endforeach()

target_link_libraries(semi_test ${ARMADILLO_LIBRARIES})

#if(OPENMP_FOUND)
#    set(TESTS_COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
#    set(TESTS_LINK_FLAGS ${OpenMP_CXX_FLAGS})
#endif(OPENMP_FOUND)

foreach(TEST ${TESTS})
    set_target_properties(${TEST} PROPERTIES
        COMPILE_FLAGS "${TESTS_COMPILE_FLAGS}"
        LINK_FLAGS "${TESTS_LINK_FLAGS}"
        BUILD_WITH_INSTALL_RPATH TRUE INSTALL_RPATH_USE_LINK_PATH TRUE)
    if(APPLE)
        set_target_properties(${TEST} PROPERTIES
            MACOSX_RPATH TRUE
            INSTALL_RPATH "${ARMADILLO_LIBRARY_DIRS}")
    endif(APPLE)
    target_link_libraries(${TEST} io libsemi ${ARMADILLO_LIBRARIES})
endforeach()
foreach(TEST ${TESTS})
    add_test(NAME ${TEST} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND ${TEST})
endforeach()

