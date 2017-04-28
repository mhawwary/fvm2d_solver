
#include"SimCase.hpp"
//#include"general_tools.h"


int main(int argc, char** argv){


    if (argc < 2) {
        FatalError_exit("ERROR: No input file specified ... ");
        return(0);
    }

    std::string input_fname;     // input file name

    input_fname = argv[argc-1];  // input file name with directory

    SimCase Case;

    Case.setup(input_fname);  // Setup and parse input parameters

    Case.InitSim();           // Preprocessing steps

    Case.RunSim();            // Main Solution Loop

    Case.PostProcess();       // Dumping Simulation Post Processing data

    return 0;
}


