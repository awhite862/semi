set(SRC
    test.cpp
)

add_library(libsemi ${SRC})
include_directories(${SEMI_DIR}/src)
