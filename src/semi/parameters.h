#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <map>
#include <sstream>

namespace Semi {

/** \brief A structure for storing calculation parameters of arbitrary type.
 **/
class parameters {
private:
    typedef std::map<std::string, std::string> map_type;
private:
    map_type m_data;

public:
    /** \brief Default constructor. **/
    parameters() { }

    /** \brief Construct from input std::map object.
        \param in_map Input map.
     **/
    parameters(const map_type &in_map) : m_data(in_map) { }

    /** \brief Returns true if the parameter structure is empty. **/
    bool is_empty() const {return m_data.empty();}

    /** \brief Set a default value if key does not exist.
        \param key Key specifying parameter to set.
        \param value Value of the defualt.
     **/
    template <typename T>
    void set_default(std::string &key, T value) {
        typename map_type::iterator i = m_data.find(key);
        if (i == m_data.end()) {
            std::stringstream ss;
            ss << value; 
            std::string svalue = ss.str();
            m_data[key] = svalue;
        }
    }

    /** \brief Add a value or overwrite existing value.
        \param key Key specifying parameter to set.
        \param value Value of parameter.
     **/
    void add(const std::string &key, const std::string &value) {m_data[key] = value;}
    
    /** \brief Get value of particular parameter
        \tparam Type of output.
        \param Key specifying parameter.
        \return Value
     **/
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


} // namesapce semi

#endif // PARAMETERS_H
