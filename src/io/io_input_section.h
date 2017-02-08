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

namespace io {

class input_section {
private:
    typedef std::map<std::string, std::string> map_type;
private:
    map_type m_data; 

public:
    input_section();
    void add(std::string &key, std::string &value);
    void add(std::string &input);
    void print(std::ostream &os);
};

} // namespace io 

#endif // IO_INPUT_SECTION_H
