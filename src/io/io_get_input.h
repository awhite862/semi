#ifndef IO_GET_INPUT_H
#define IO_GET_INPUT_H

#include "io_input_section.h"
#include "io_molecule.h"

namespace io {

void get_input_sections(
    std::istream &is, 
    std::map<std::string, molecular_input> &min,
    std::map<std::string, input_section> &input);

void print_input_sections(
    std::ostream &os, 
    std::map<std::string, molecular_input> &min,
    std::map<std::string, input_section> &input);

} // namespace io

#endif // IO_GET_INPUT_H
