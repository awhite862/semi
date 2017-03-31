#include "io_input_section.h"
#include <iomanip>

namespace io {

input_section::input_section() {  } 
    
const std::map<std::string, std::string>& input_section::get() const {
    return m_data;
}
    
bool input_section::is_empty() const {return m_data.empty();}

void input_section::add(std::string &key, std::string &value) {
    m_data[key] = value;
}

void input_section::print(std::ostream &os) const {
    for (typename map_type::const_iterator ii = m_data.begin();
        ii != m_data.end(); ++ii) {

        os << " " << ii->first << " " << ii->second << std::endl;
    }
}

} // namespace io
