#include <fstream>
#include <cstdlib> // for EXIT_FAILURE
#include <io/io_get_input.h>
#include <semi/semi_api.h>

//TODO move this to utilities
double get_charge(std::string s) {
    if (isdigit(s[0])) {
        double charge;
        std::stringstream ss;
        ss << s;
        ss >> charge;
        return charge;
    }
    else {
        if (s == "h" ) return 1.0;
        else if (s == "he") return 2.0;
        else if (s == "li") return 3.0;
        else if (s == "be") return 4.0;
        else if (s == "b" ) return 5.0;
        else if (s == "c" ) return 6.0;
        else if (s == "n" ) return 7.0;
        else if (s == "o" ) return 8.0;
        else if (s == "f" ) return 9.0;
        else if (s == "ne") return 10.0;
        else return 0.0;
    }
}

void build_input(
    const std::map<std::string, io::molecular_input> &min,
    const std::map<std::string, io::input_section> &inputs,
    Semi::input &in) {

    typedef std::map<std::string, io::input_section> imap_type;
    typedef std::map<std::string, io::molecular_input> mmap_type;

    size_t nmol = min.size();
    if (nmol == 1) {
        typename mmap_type::const_iterator ii = min.find("molecule");
        if (ii == min.end()) {
            throw std::runtime_error("Could not find 'molecule' section");
        }
        std::vector<Semi::Atom> avec;
        for (size_t i = 0; i < ii->second.size(); i++) {
            double charge = get_charge(ii->second.gets(i));
            Semi::Atom A(
                ii->second.getx(i), 
                ii->second.gety(i), 
                ii->second.getz(i), 
                charge);

            avec.push_back(A);
        }
        in.mol = Semi::Molecule(avec);
    }
    else {
        std::stringstream ss;
        ss << "Incorrect molecule specification: " << nmol << " molecules were specified.";
        throw std::runtime_error(ss.str()); 
    }
    
    // parse control section
    {
        typename imap_type::const_iterator ii = inputs.find("control");
        if (ii == inputs.end()) {
            throw std::runtime_error("Could not find 'control' section");
        }
        const io::input_section &cont = ii->second;
        std::string method = cont.get_value<std::string>("method");
        if (method == "huckel") {
            in.ctype = Semi::HUCKEL;
        }
        else if (method == "cndo") {
            in.ctype = Semi::CNDO;
        }
        else {
            throw std::runtime_error("Unrecognized method: " + method);
        }
    }

    // get huckel section if applicable
    if (inputs.count("huckel") != 0) {
        typename imap_type::const_iterator ih = inputs.find("huckel");
        in.huckel_params = Semi::parameters(ih->second.get());
    }

    // get CNDO section if applicable
    if (inputs.count("cndo") != 0) {
        typename imap_type::const_iterator ic = inputs.find("cndo");
        in.huckel_params = Semi::parameters(ic->second.get());
    }
}

int main(int argc, char * argv[]) {

    try {
        if (argc < 2) {
            std::cout << "No input file was specified!" << std::endl;
            return EXIT_FAILURE;
        }
        else if (argc < 4) {

            // Parse input file
            std::ifstream fin(argv[1]);
            std::map<std::string, io::molecular_input> min;
            std::map<std::string, io::input_section> inputs;
            io::get_input_sections(fin, min, inputs);

            // Build input for Semi
            Semi::input in;
            build_input(min, inputs, in);
            
            // Run calculation printing output to std::cout
            if (argc == 2) {
                io::print_input_sections(std::cout, min, inputs);
                Semi::output out;
                try {
                    Semi::run_semi(in, out, std::cout);
                }
                catch (std::exception &e) {
                    std::cout << "Exception in Semi: " << e.what() << std::endl;
                    return EXIT_FAILURE;
                }
            }

            // Run calculation printing output to file
            else if (argc == 3) {
                std::ofstream fout(argv[2]);
                io::print_input_sections(fout, min, inputs);
                Semi::output out;
                try {
                    Semi::run_semi(in, out, fout);
                }
                catch (std::exception &e) {
                    fout << "Exception in Semi: " << e.what() << std::endl;
                    return EXIT_FAILURE;
                }
            }
        }
        else {
            std::cout << "Too many input parameters." << std::endl;
            return EXIT_FAILURE;
        }
    }
    catch (std::exception &e) {
        std::cout << "Exception: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    catch (...) {
        std::cout << "Unrecognized exception!!!" << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
