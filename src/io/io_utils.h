#ifndef IO_UTILS_H
#define IO_UTILS_H

#include <string>

namespace io {

/** \brief Clean up input string
    
    In particular,
        -- remove leading and trailing whitespace
        -- remove all characters following a comment '#' symbol
    \param[in,out] str String to be cleaned.
 **/
void clean_string(std::string &str);

void replace_tabs(std::string &str);

void trim_whitespace(std::string &str);

void trim_comment(std::string &str, const std::string &com = "#"); 

} // namespace io 

#endif // IO_UTILS
