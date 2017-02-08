//
//      Copyright Alec White 2017
// Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE.txt or copy at
//      http://www.boost.org/LICENSE_1_0.txt)
//
#ifndef IO_INPUT_SECTION_H
#define IO_INPUT_SECTION_H

#include <map>
#include <string>
#include <iostream>
#include <sstream>

namespace io {

class input_section {
private:
    typedef std::map<std::string, std::string> map_type;
private:
    map_type m_data; 

public:
    input_section();
    void add(std::string &key, std::string &value);
    void print(std::ostream &os) const;
    bool is_empty() const;
    
    template<typename T>
    T get_value(const std::string &key) {
        if(m_data.count(key) == 0) {
            throw std::logic_error(
                "input_section::get_value(): parameter not found " + key);
        }
        const std::string val = m_data[key];
        T out;
        std::stringstream ss;
        ss << val; ss >> out;
        return out;
    }
};

void get_input_sections(
    std::istream &is, 
    std::map<std::string, input_section> &input);

void print_input_sections(
    std::ostream &os, 
    std::map<std::string, input_section> &input);

} // namespace io 

#endif // IO_INPUT_SECTION_H
