#include "io_input_section.h"
#include <iomanip>

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


} // namespace io
