/*
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

Function to set up data directories and copy input files into them to retain
a record of which settings and data were used for the run. If a directory
already exists with the same name as the output directory being set up, the old
one is renamed to avoid accidental deletion. Any previous existing version of
this 'previous output' directory is deleted. */

/////////////////////////////////////////////////////
/*
Copyright (C) 2023, Daniel Duffy, daniellouisduffy@gmail.com. All rights reserved.
Please cite Daniel Duffy and John S. Biggins if you 
use any part of this code in work that you publish or distribute.

This file is part of MorphoShell.

MorphoShell is distributed under the terms of the Cambridge Academic
Software License (CASL). You should have received a copy of the license
along with MorphoShell. If not, contact Daniel Duffy, daniellouisduffy@gmail.com.
*/
/////////////////////////////////////////////////////

// Turn Eigen bounds checking off for speed.
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <filesystem> // Note this requires g++ version >=8

#include "extract_filename.hpp"
#include "initialize_directories.hpp"
#include "Out_Stream_Class.hpp"
#include "Stuff_Class.hpp"


void initialize_directories(
    Stuff_Class &stuff,
    std::string &first_log_entry,
    Out_Stream_Class &log){

    // Get full path of current working directory.
    std::filesystem::path cwd = std::filesystem::current_path();
    first_log_entry += "Current working directory: " + cwd.string() + "\n";

    // Create strings containing names for the output directory and a potential
    // backup for the previous one.
    std::string output_dir = "output_using_" + extract_filename(stuff.settings_filepath);
    std::string prev_output_dir = "prev_" + output_dir;

    // If previous version of outputfiles_using... directory already exists for the
    // current settings file name, save (rename) it as prev_output...
    // (**OVERWRITTEN EVERY TIME THIS IS DONE!**).
    std::filesystem::path to_rename(cwd / output_dir);

    if( std::filesystem::exists(to_rename) ){

        // Delete any existing prev_output... directory corresponding to
        // the current settings file (i.e. if the code saved you from accidentally
        // overwriting output_using_... the first time, but you forgot to rename
        // the prev_output... directory that saved you: unlucky!).
        std::filesystem::path to_delete(cwd / prev_output_dir);
        if( std::filesystem::exists(to_delete) ){
            std::filesystem::remove_all(to_delete);
            first_log_entry +=
            "Deleted existing output directory from run-before-last: " + prev_output_dir + "\n";
        }

        // Do renaming of output_using...
        std::filesystem::rename(to_rename, cwd / prev_output_dir);
        first_log_entry += "Found a previous output directory for these settings.\n";
        first_log_entry +=
        "The previous version was saved this time, but you may want to move it "
        "out of the working directory, \n because it will be overwritten next time "
        "the code is run with this settings file.\n";
    }
    else{
        first_log_entry += "No previous output directory was found that needed saving.\n";
    }


    // Create new directory to hold output .vtk files, and copy settings and data
    // files into a subdirectory of that directory.
    bool mk_output_dir =  std::filesystem::create_directory(std::filesystem::path(cwd / output_dir));
    std::filesystem::path input_files_dir(cwd / (output_dir+"/input_files_used"));
    bool mk_input_files_dir =  std::filesystem::create_directory(input_files_dir);

    bool cp_settings = std::filesystem::copy_file(cwd / stuff.settings_filepath, input_files_dir / extract_filename(stuff.settings_filepath));
    bool cp_ref_data = std::filesystem::copy_file(cwd / stuff.ref_data_filepath, input_files_dir / extract_filename(stuff.ref_data_filepath));
    bool cp_ansatz = true; // Default true so no error thrown if no ansatz given.
    if( stuff.ansatz_file_was_given ){
        cp_ansatz = std::filesystem::copy_file(cwd / stuff.ansatz_filepath, input_files_dir / extract_filename(stuff.ansatz_filepath));
    }

    if(mk_output_dir == false || mk_input_files_dir == false || cp_settings == false || cp_ref_data == false || cp_ansatz == false){
        // This exception on top of the std::filesystem ones may be overkill.
        throw std::runtime_error("Error: setting up output directory failed");
    }



    // We also set up the log.txt file, and give it its first entry.
    log.set_output_filepath(output_dir + "/log.txt");
    // Check we can create a log file with the chosen name.
    try{
        log.open();
        log.close();
    }
    catch(const std::runtime_error &error){
        std::cerr << error.what() << std::endl;
        }
    log << std::defaultfloat << std::setprecision(6);
    log << first_log_entry << std::endl;


    // We also set up an auxiliary data file, which is handy 
    // for outputting 'extra' things, e.g. the dial_factor.
    stuff.aux_output_filepath = output_dir + "/aux_output.txt";
    std::ofstream aux_output;
    aux_output.open(stuff.aux_output_filepath, std::ofstream::app);
    if( !aux_output ){
        throw std::runtime_error("Error: Problem creating aux_output filestream.");
    }
    aux_output.close();


    // We store the path to the output dir, 
    // converted to a std::string.
    stuff.output_dir_path = std::string(cwd / output_dir);
}
