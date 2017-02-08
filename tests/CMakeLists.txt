set(TESTS 
    arma_test
    io_utils_test
    framework_test
)
foreach(TEST ${TESTS})
    add_executable(${TEST} ${TEST}.cpp)
endforeach()
include_directories(${SEMI_DIR}/src)

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
    target_link_libraries(${TEST} io
        ${ARMADILLO_LIBRARIES})
endforeach()
foreach(TEST ${TESTS})
    add_test(NAME ${TEST} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND ${TEST})
endforeach()

