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
    const map_type& get() const;
    
    template<typename T>
    T get_value(const std::string &key) const {
        if(m_data.count(key) == 0) {
            throw std::logic_error(
                "input_section::get_value(): parameter not found " + key);
        }
        const std::string val = m_data.find(key)->second;
        T out;
        std::stringstream ss;
        ss << val; ss >> out;
        return out;
    }
};

} // namespace io 

#endif // IO_INPUT_SECTION_H
