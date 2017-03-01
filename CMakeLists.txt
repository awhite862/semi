cmake_minimum_required(VERSION 2.8.4)
project(SEMI CXX)

set(SEMI_DIR "${SEMI_SOURCE_DIR}")

include(CheckCXXSourceCompiles)
include(CheckIncludeFileCXX)
include("${SEMI_DIR}/cmake/SEMIDeps.txt")

#find_package(Doxygen)

include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${SEMI_DIR}/src)
set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}")

enable_testing()
add_subdirectory(src/io)
add_subdirectory(src/semi)
add_subdirectory(tests)
add_executable(semi++ ${SEMI_DIR}/src/main.cpp)
target_link_libraries(semi++ io libsemi)

#if(DOXYGEN_FOUND AND (DOXYGEN_VERSION VERSION_GREATER "1.8"))
#    configure_file(
#        "${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in"
#        "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile" @ONLY)
#    add_custom_target(docs
#        "${DOXYGEN_EXECUTABLE}" "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile"
#        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
#        COMMENT "Doxygen documentation" VERBATIM)
#endif()
