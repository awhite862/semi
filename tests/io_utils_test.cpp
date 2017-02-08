#include <iostream>
#include <io/io_utils.h>

using namespace io;

int run_clean_test() {
    std::string t = "  this is\ta test   ";
    const std::string o = "this is a test";
    clean_string(t);
    if (t != o) return 1;
    return 0;
}

int run_comment_test() {
    std::string t = "  this is a test # This is a comment   ";
    const std::string o = "this is a test";
    clean_string(t);
    if (t != o) {
        std::cout << " expected: " << o << std::endl;
        std::cout << " found: " << t << std::endl;
        return 1;
    }
    return 0;
}


int main() {
    return 

    run_clean_test() |
    run_comment_test() |

    0;
}
