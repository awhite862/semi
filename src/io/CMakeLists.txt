set(SRC
    io_get_input.cpp
    io_input_section.cpp
    io_molecule.cpp
    io_utils.cpp
)

add_library(io ${SRC})
include_directories(${SEMI_DIR}/src)
