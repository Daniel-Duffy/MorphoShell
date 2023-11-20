/*
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

Function to read in settings file and put these
these settings in a class which is then returned.
*/

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

#include <iostream> // Used for outputting messages to the terminal
#include <fstream> // Used for file I/O.
#include <stdexcept> // Used for exception (error) handling
#include <string> // Used for creating directory for output data
#include <sstream>
#include <vector>
#include <list>
#include <cstddef> // Defines std::size_t
#include <algorithm> // For std::replace and std::remove_if etc
#include <cctype> // For std::isspace

#include "read_settings.hpp"
#include "Stuff_Class.hpp"
#include "Out_Stream_Class.hpp"

std::vector<std::size_t> substr_positions_in_str(const std::string &str, const std::string &sub);
bool boolstring_to_bool(const std::string &boolstring);

void read_settings(
    Stuff_Class &stuff,
    Out_Stream_Class &log){

std::string file_contents = "";

std::vector<std::string> settings_vec;
std::vector<std::string> setting_names_vec;


std::ifstream settings_file(stuff.settings_filepath);
if( !settings_file ){
    throw std::runtime_error("Error: Problem opening settings file.");
}

std::stringstream temp;
temp << settings_file.rdbuf();
file_contents = temp.str();
// Get rid of any carriage returns, which I think are basically a Windows thing.
std::replace(file_contents.begin(), file_contents.end(), '\r', '\n');
file_contents.append("\n"); // This makes the next step a little simpler.

// If any "//" is encountered, erase from there to the end of that line.
// In this way we support commenting out things in the settings file with //.
while( true ){

    long unsigned int first_comment_pos = file_contents.find("//");

    // Break if no comment was found.
    if( first_comment_pos == std::string::npos ){
        break;
    }

    long unsigned int next_newline_pos = file_contents.find("\n", first_comment_pos);

    // Erase the comment.
    file_contents.erase(first_comment_pos, (next_newline_pos - first_comment_pos) + 1);
}


// Now remove all white space. After doing this file_contents
// should just look like setting1=1.9;setting2=7;setting3=hello;
// for example.
// Here we are using the 'erase-remove idiom' --- see wiki 
// or std::remove_if documentation.
file_contents.erase(std::remove_if(file_contents.begin(), file_contents.end(), [](unsigned char x) { return std::isspace(x); }), file_contents.end());
//std::cout << file_contents << std::endl;

// In the remaining string, find positions of all equals signs and semicolons.
std::vector<std::size_t> equalsign_positions = substr_positions_in_str("=", file_contents);
std::vector<std::size_t> semicolon_positions = substr_positions_in_str(";", file_contents);

if( equalsign_positions.size() != semicolon_positions.size() ){
    throw std::runtime_error("Error: After removing all comments, the settings file did not have equal numbers of = and ; \n The most common mistake is to forget a semicolon.");
}
std::vector<std::size_t> vec_that_should_be_ascending;
for(std::size_t j = 0; j < equalsign_positions.size(); ++j){
    vec_that_should_be_ascending.push_back(equalsign_positions.at(j));
    vec_that_should_be_ascending.push_back(semicolon_positions.at(j));
}
if( not std::is_sorted(vec_that_should_be_ascending.begin(), vec_that_should_be_ascending.end()) ){
    throw std::runtime_error("Error: After removing all comments, the = and ; in the settings file did not occur in the correct order.");
}

size_t cursor = 0;
size_t k = 0;
while( true ){
    setting_names_vec.push_back(file_contents.substr( cursor, equalsign_positions.at(k) - cursor ));
    cursor = equalsign_positions.at(k) + 1;
    settings_vec.push_back(file_contents.substr( cursor, semicolon_positions.at(k) - cursor ));
    cursor = semicolon_positions.at(k) + 1;
    ++k;
    if( k == equalsign_positions.size() ){
        break;
    }
}

// Print std::vector (stackoverflow 10750057).
//for( const auto& v: setting_names_vec ){
//    std::cout << v << std::endl;
//}
//for( const auto& v: settings_vec ){
//    std::cout << v << std::endl;
//}

// Now transfer the settings to the stuff_struct.

// If you add a setting, you must add a corresponding line 
// at the correct point in the following sequence (and 
// the order of the settings in the file matters!).
int j = -1;
//
++j;
if(setting_names_vec.at(j) != "outputs_per_char_long_time"){throw std::runtime_error("Error: something went wrong when reading the outputs_per_char_long_time setting.");}
stuff.outputs_per_char_long_time = std::stod(settings_vec.at(j));
//
//++j;
//if(setting_names_vec.at(j) != "inty"){throw std::runtime_error("Error: something went wrong when reading the inty setting.");}
//stuff.inty = std::stoi(settings_vec.at(j));
//
//++j;
//if(setting_names_vec.at(j) != "stringy"){throw std::runtime_error("Error: something went wrong when reading the stringy setting.");}
//stuff.stringy = settings_vec.at(j);
//
++j;
if(setting_names_vec.at(j) != "ref_thickness_if_uniform"){throw std::runtime_error("Error: something went wrong when reading the ref_thickness_if_uniform setting.");}
stuff.ref_thickness_if_uniform = std::stod(settings_vec.at(j));
//
++j;
if(setting_names_vec.at(j) != "timestep_prefactor"){throw std::runtime_error("Error: something went wrong when reading the timestep_prefactor setting.");}
stuff.timestep_prefactor = std::stod(settings_vec.at(j));
//
++j;
if(setting_names_vec.at(j) != "dialling_from_ansatz_rather_than_ref_state"){throw std::runtime_error("Error: something went wrong when reading the dialling_from_ansatz_rather_than_ref_state setting.");}
stuff.dialling_from_ansatz_rather_than_ref_state = boolstring_to_bool(settings_vec.at(j));
//
++j;
if(setting_names_vec.at(j) != "dial_factors_to_hold_at_spacing"){throw std::runtime_error("Error: something went wrong when reading the dial_factors_to_hold_at_spacing setting.");}
stuff.dial_factors_to_hold_at_spacing = std::stod(settings_vec.at(j));
//
++j;
if(setting_names_vec.at(j) != "dial_phase_time_prefactor"){throw std::runtime_error("Error: something went wrong when reading the dial_phase_time_prefactor setting.");}
stuff.dial_phase_time_prefactor = std::stod(settings_vec.at(j));
//
++j;
if(setting_names_vec.at(j) != "time_between_equil_checks_prefactor"){throw std::runtime_error("Error: something went wrong when reading the time_between_equil_checks_prefactor setting.");}
stuff.time_between_equil_checks_prefactor = std::stod(settings_vec.at(j));
//
++j;
if(setting_names_vec.at(j) != "damping_prefactor_1"){throw std::runtime_error("Error: something went wrong when reading the damping_prefactor_1 setting.");}
stuff.damping_prefactor_1 = std::stod(settings_vec.at(j));
//
++j;
if(setting_names_vec.at(j) != "damping_prefactor_2"){throw std::runtime_error("Error: something went wrong when reading the damping_prefactor_2 setting.");}
stuff.damping_prefactor_2 = std::stod(settings_vec.at(j));
//
++j;
if(setting_names_vec.at(j) != "num_steps_to_use_gradient_descent_dynamics_for_before_switching_to_newtonian_dynamics"){throw std::runtime_error("Error: something went wrong when reading the num_steps_to_use_gradient_descent_dynamics_for_before_switching_to_newtonian_dynamics setting.");}
stuff.num_steps_to_use_gradient_descent_dynamics_for_before_switching_to_newtonian_dynamics = std::stoi(settings_vec.at(j));
//
++j;
if(setting_names_vec.at(j) != "force_to_char_force_ratio_equil_threshold"){throw std::runtime_error("Error: something went wrong when reading the force_to_char_force_ratio_equil_threshold setting.");}
stuff.force_to_char_force_ratio_equil_threshold = std::stod(settings_vec.at(j));
//
++j;
if(setting_names_vec.at(j) != "speed_to_char_speed_ratio_equil_threshold"){throw std::runtime_error("Error: something went wrong when reading the speed_to_char_speed_ratio_equil_threshold setting.");}
stuff.speed_to_char_speed_ratio_equil_threshold = std::stod(settings_vec.at(j));
//
++j;
if(setting_names_vec.at(j) != "def_poisson_ratio"){throw std::runtime_error("Error: something went wrong when reading the def_poisson_ratio setting.");}
stuff.def_poisson_ratio = std::stod(settings_vec.at(j));
//
++j;
if(setting_names_vec.at(j) != "ref_shear_modulus_if_uniform"){throw std::runtime_error("Error: something went wrong when reading the ref_shear_modulus_if_uniform setting.");}
stuff.ref_shear_modulus_if_uniform = std::stod(settings_vec.at(j));
//
++j;
if(setting_names_vec.at(j) != "ref_density"){throw std::runtime_error("Error: something went wrong when reading the ref_density setting.");}
stuff.ref_density = std::stod(settings_vec.at(j));
//
++j;
if(setting_names_vec.at(j) != "bend_stiffening_scale_factor"){throw std::runtime_error("Error: something went wrong when reading the bend_stiffening_scale_factor setting.");}
stuff.bend_stiffening_scale_factor = std::stod(settings_vec.at(j));
//
++j;
if(setting_names_vec.at(j) != "patch_mat_dimless_conditioning_thresh"){throw std::runtime_error("Error: something went wrong when reading the patch_mat_dimless_conditioning_thresh setting.");}
stuff.patch_mat_dimless_conditioning_thresh = std::stod(settings_vec.at(j));
//
++j;
if(setting_names_vec.at(j) != "slide_stiffness_prefactor"){throw std::runtime_error("Error: something went wrong when reading the slide_stiffness_prefactor setting.");}
stuff.slide_stiffness_prefactor = std::stod(settings_vec.at(j));
//
++j;
if(setting_names_vec.at(j) != "slide_friction_coefficent"){throw std::runtime_error("Error: something went wrong when reading the slide_friction_coefficent setting.");}
stuff.slide_friction_coefficent = std::stod(settings_vec.at(j));
//
++j;
if(setting_names_vec.at(j) != "slide_speed_prefactor"){throw std::runtime_error("Error: something went wrong when reading the slide_speed_prefactor setting.");}
stuff.slide_speed_prefactor = std::stod(settings_vec.at(j));
//
++j;
if(setting_names_vec.at(j) != "using_initial_glass_slide_z_coords_from_settings_file"){throw std::runtime_error("Error: something went wrong when reading the using_initial_glass_slide_z_coords_from_settings_file setting.");}
stuff.using_initial_glass_slide_z_coords_from_settings_file = boolstring_to_bool(settings_vec.at(j));
//
++j;
if(setting_names_vec.at(j) != "initial_upper_slide_z_coord"){throw std::runtime_error("Error: something went wrong when reading the initial_upper_slide_z_coord setting.");}
stuff.initial_upper_slide_z_coord = std::stod(settings_vec.at(j));
//
++j;
if(setting_names_vec.at(j) != "initial_lower_slide_z_coord"){throw std::runtime_error("Error: something went wrong when reading the initial_lower_slide_z_coord setting.");}
stuff.initial_lower_slide_z_coord = std::stod(settings_vec.at(j));
//



if( stuff.dialling_from_ansatz_rather_than_ref_state ){
    if( stuff.ansatz_file_was_given ){
        log << "You have chosen to dial abar and bbar from the a and b of the ansatz. "
                "This will override some things: dial_factor_to_start_from will be set "
                "to 0, and a and b will be dialled linearly from the ansatz a and b to "
                "the final (dial_factor=1) abar and bbar." << std::endl;
    }
    else{
        throw std::runtime_error("Error: You cannot use the dialling_from_ansatz_rather_than_ref_state setting if you have not supplied an ansatz file.");
        }
}

if( stuff.dial_factors_to_hold_at_spacing > 1.0 ){
    throw std::runtime_error("Error: The largest permitted value of dial_factors_to_hold_at_spacing is 1.0, because dial_factor never goes above 1.0, so anything larger would be meaningless.");
}

if( stuff.num_steps_to_use_gradient_descent_dynamics_for_before_switching_to_newtonian_dynamics > 0){
    stuff.using_gradient_descent_dynamics = true;
    log << "You have chosen to use gradient descent dynamics for " << stuff.num_steps_to_use_gradient_descent_dynamics_for_before_switching_to_newtonian_dynamics << " steps before switching to Newtonian dynamics." << std::endl;
}
else{
    stuff.using_gradient_descent_dynamics = false;
}


}


// Function to take a std:string str, and find all the positions of
// (first characters of) the occurences of a substring sub in str.
std::vector<std::size_t> substr_positions_in_str(const std::string &sub, const std::string &str){

    std::vector<std::size_t> positions;
    std::size_t pos = str.find(sub, 0);

    while(pos != std::string::npos){

    positions.push_back(pos);
    pos = str.find(sub,pos+1);
    }

    return positions;
}



// Function to take a std::string "true" and convert it to
// an actual bool, and similarly for false.
bool boolstring_to_bool(const std::string &boolstring){
    if( boolstring == "true" ){
        return true;
    }
    else if( boolstring == "false" ){
        return false;
    }
    else{
        throw std::runtime_error("Error: you passed a std::string to boolstring_to_bool that read neither true nor false.");
    }
}
