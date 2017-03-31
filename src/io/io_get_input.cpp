#include <algorithm>
#include "io_get_input.h"
#include "io_utils.h"

namespace io {

namespace {

void expand_molecular_abbreviations(std::string &str) {
    if (str == "mol") str = "molecule";
}

} // unnamed namespace 

void get_input_sections(
    std::istream &is, 
    std::map<std::string, molecular_input> &min,
    std::map<std::string, input_section> &input) {

    std::string name;
    while ( std::getline(is, name) ) {
        clean_string(name);
        if (!name.empty() && name[0] == '$') {
            name = name.substr(1,name.length());
            clean_string(name);
            std::transform(name.begin(), name.end(), name.begin(), ::tolower);
            input[name] = input_section();
            std::string temp;
            while ( std::getline(is, temp) ) {
                clean_string(temp);
                std::transform(temp.begin(), temp.end(), temp.begin(), ::tolower);
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
        else if (!name.empty() && name[0] == '@') {
            name = name.substr(1,name.length());
            clean_string(name);
            std::transform(name.begin(), name.end(), name.begin(), ::tolower);
            expand_molecular_abbreviations(name);
            min[name] = molecular_input();
            std::string temp;
            while ( std::getline(is, temp) ) {
                clean_string(temp);
                std::transform(temp.begin(), temp.end(), temp.begin(), ::tolower);
                if (!temp.empty() && temp[0] == '@') break;
                min[name].add(temp);
            }
        }
    }
}

void print_input_sections(
    std::ostream &os, 
    std::map<std::string, molecular_input> &min,
    std::map<std::string, input_section> &input) {

    typedef std::map<std::string, molecular_input>::iterator imtype;
    for (imtype ii = min.begin(); ii != min.end(); ++ii) {
        os << ii->first << ": " << std::endl;
        ii->second.print(os);
    }

    typedef std::map<std::string, input_section>::iterator istype;
    for (istype ii = input.begin(); ii != input.end(); ++ii) {
        os << ii->first << ": " << std::endl;
        ii->second.print(os);
    }
}

} // namespace io
