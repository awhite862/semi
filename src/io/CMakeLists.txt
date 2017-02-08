set(SRC
    io_input_section.cpp
    io_utils.cpp
)

add_library(io ${SRC})
include_directories(${SEMI_DIR}/src)
