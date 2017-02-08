//
//      Copyright Alec White 2017
// Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE.txt or copy at
//      http://www.boost.org/LICENSE_1_0.txt)
//
#include "io_utils.h"

namespace io {

void clean_string(std::string &str) {
    trim_comment(str);
    replace_tabs(str);
    trim_whitespace(str);
}

// TODO replace all this with regex
void replace_tabs(std::string &str) {
    const std::string ts = "\t";
    const std::string ss = " ";
    const size_t max = 2000;
    size_t count = 0;

    size_t found = str.find(ts);
    while (found != std::string::npos && count < max) {
        str.replace(found, ts.length(), ss);
        found = str.find(ts);
        count++;
    }
}

void trim_whitespace(std::string &str) {
    const char w = ' ';
    size_t ifirst = str.find_first_not_of(w);
    ifirst = (ifirst == std::string::npos) ? str.length() : ifirst;
    size_t ilast = str.find_last_not_of(w);
    str = str.substr(ifirst, ilast - ifirst + 1);
}

void trim_comment(std::string &str, const std::string &com) {
    const size_t ifirst = 0;
    size_t ilast = str.find(com);
    ilast = (ilast == std::string::npos) ? str.length() : ilast;
    str = str.substr(ifirst, ilast - ifirst);
}

} // namespace io 
