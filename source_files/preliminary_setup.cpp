/*
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

Function to set up data directories and copy input files into them to retain
a record of which settings and data were used for the run. If a directory
already exists with the same name as the output directory being set up, the old
one is renamed to avoid accidental deletion. Any previous existing version of
this 'previous output' directory is deleted. */

#include <iostream>
#include <string>
#include <filesystem> // Note this requires g++ version >=8

#include "extract_filename.hpp"
#include "initialize_directories.hpp"
#include "Out_Stream_Class.hpp"
#include "Stuff_Class.hpp"
#include "get_real_time.hpp"

#include "preliminary_setup.hpp"



void preliminary_setup(
    int argc, 
    char *argv[],
    Stuff_Class &stuff,
    Out_Stream_Class &log
 ){



// Check command line arguments, with ansatz file optional, 
// and store the file paths.
if( (argc != 3) && (argc != 4) )
    {std::cerr << "Error: Please provide at least one settings file path "
    "followed by one data file path as command line arguments.\n"
    "A path to a file with ansatz degrees of freedom can be given as a third "
    "argument if desired." << std::endl;
    }
stuff.settings_filepath = std::string(argv[1]);
stuff.ref_data_filepath = std::string(argv[2]);
if( argc == 3 ){
    stuff.ansatz_filepath = "no_ansatz_file";
    stuff.ansatz_file_was_given = false;
}
else{
    stuff.ansatz_filepath = std::string(argv[3]);
    stuff.ansatz_file_was_given = true;
}

// The log file will store a copy of the std::cout and std::cerr
// output. It's stream will be passed by reference to any functions that 
// want to output to the log. The stream must be opened before writing to it, 
// and it's good practice to close it immediately afterwards.
// Store some handy things that will go in the log.
std::string first_log_entry;
first_log_entry += "Simulation run began at: " + get_real_time() + "\n";
first_log_entry += "The command that was run was\n";
for(int i = 0; i < argc; ++i){
    first_log_entry +=  argv[i];
    first_log_entry +=  " ";
}
first_log_entry += "\n";
first_log_entry += "not including any openMP part it might have had in front,\n";
first_log_entry += "which is usually something like\n";
first_log_entry += "OMP_PLACES=cores OMP_PROC_BIND=spread OMP_NUM_THREADS=5\n";
first_log_entry += "\n\n";

// Now set up output directory, log file etc.
try{ 
    initialize_directories( 
    stuff,
    first_log_entry,
    log ); }
catch(const std::runtime_error &error){
    std::cerr << error.what() << std::endl;
    std::cerr << "Error: initialize_directories failed." << std::endl;
    }



}
