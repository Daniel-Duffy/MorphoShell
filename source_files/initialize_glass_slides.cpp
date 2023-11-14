/*
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

*/

//Turn Eigen bounds checking off for speed.
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <vector>
#include <Eigen/Dense>
#include <string>

#include "initialize_glass_slides.hpp"
#include "Stuff_Class.hpp"
#include "find_extremum_of_doubles.hpp"


void initialize_glass_slides(
    const Eigen::Matrix<double,Eigen::Dynamic,3> &node_positions,
    Stuff_Class &stuff){

    // Can specify initial slide z coords directly in settings file.
    if( stuff.using_initial_glass_slide_z_coords_from_settings_file ){
    stuff.upper_slide_z_coord = stuff.initial_upper_slide_z_coord;
    stuff.lower_slide_z_coord = stuff.initial_lower_slide_z_coord;
    }
    // Otherwise set the slides to be aligned with the top and bottom 
    // of the shell.
    else{
        std::vector<double> node_z_values;
        for(int n = 0; n < stuff.num_nodes; ++n){
            node_z_values.push_back(node_positions(n,2));
        }
        stuff.upper_slide_z_coord = min_or_max_of_doubles(node_z_values, std::string("max"));
        stuff.lower_slide_z_coord = min_or_max_of_doubles(node_z_values, std::string("min"));
    }    
}
