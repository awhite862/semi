#include <fstream>
#include <io/io_input_section.h>

int main(int argc, char * argv[]) {
     
    try {

        if (argc < 2) {
            std::cout << "No input file was specified!" << std::endl;
            return EXIT_FAILURE;
        }
        else if (argc < 4) {
            std::ifstream fin(argv[1]);
            std::map<std::string, io::input_section> inputs;
            io::get_input_sections(fin, inputs);
            if (argc == 2) {
                io::print_input_sections(std::cout, inputs);

            }
            else if (argc == 3) {
                std::ofstream fout(argv[2]);
                io::print_input_sections(fout, inputs);
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
