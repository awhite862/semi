//
//      Copyright Alec White 2017
// Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE.txt or copy at
//      http://www.boost.org/LICENSE_1_0.txt)
//
#include "io_input_section.h"
#include <iomanip>
#include <sstream>
#include "io_utils.h"

namespace io {

input_section::input_section() {  } 

void input_section::add(std::string &key, std::string &value) {
    m_data[key] = value;
}

void input_section::add(std::string &input) {
    clean_string(input);
    std::stringstream ss;
    std::string key, value;
    ss << input;
    ss >> key; ss >> value;
    if (!key.empty() && !value.empty()) {
        m_data[key] = value;
    }
}

void input_section::print(std::ostream &os) {
    for (typename map_type::iterator ii = m_data.begin();
        ii != m_data.end(); ++ii) {

        os << ii->first << " " << ii->second << std::endl;
    }
}

} // namespace io
