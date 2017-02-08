//
//      Copyright Alec White 2017
// Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE.txt or copy at
//      http://www.boost.org/LICENSE_1_0.txt)
//
#include "io_input_section.h"
#include <iomanip>
#include "io_utils.h"

namespace io {

input_section::input_section() {  } 

void input_section::add(std::string &key, std::string &value) {
    m_data[key] = value;
}

void input_section::print(std::ostream &os) const {
    for (typename map_type::const_iterator ii = m_data.begin();
        ii != m_data.end(); ++ii) {

        os << " " << ii->first << " " << ii->second << std::endl;
    }
}
    
bool input_section::is_empty() const {return m_data.empty();}

void get_input_sections(
    std::istream &is, 
    std::map<std::string, input_section> &input) {

    std::string name;
    while ( std::getline(is, name) ) {
        clean_string(name);
        if (!name.empty() && name[0] == '$') {
            name = name.substr(1,name.length());
            input[name] = input_section();
            std::string temp;
            while ( std::getline(is, temp) ) {
                clean_string(temp);
                if (!temp.empty() && temp[0] == '$') break;
                std::string key, value;
                std::stringstream ss;
                if (!temp.empty()) {
                    ss << temp;
                    ss >> key; ss >> value;
                    input[name].add(key, value);
                }
            }
        }
    }
}

void print_input_sections(
    std::ostream &os, 
    std::map<std::string, input_section> &input) {

    typedef std::map<std::string, input_section>::iterator itype;
    for (itype ii = input.begin(); ii != input.end(); ++ii) {
        os << ii->first << ": " << std::endl;
        ii->second.print(os);
    }
}

} // namespace io
