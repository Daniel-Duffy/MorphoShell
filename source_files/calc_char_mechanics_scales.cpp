/*
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

Calculate and store characteristic force, energy, and energy density scales.

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

#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>

#include "find_extremum_of_doubles.hpp"
#include "kahan_sum.hpp"
#include "Tri_Class.hpp"
#include "Stuff_Class.hpp"
#include "Out_Stream_Class.hpp"
#include "calc_char_mechanics_scales.hpp"

void calc_char_mechanics_scales(
    const std::vector<Tri_Class> &triangles,
    Stuff_Class &stuff){

    // Calc and store minimum reference thickness.
    std::vector <double> ref_thicknesses;
    for( auto& tri: triangles ){
        ref_thicknesses.push_back(tri.ref_thickness);
    }
    stuff.min_ref_thickness = min_or_max_of_doubles(ref_thicknesses, std::string("min"));

    // Calc and store minimum reference shear_modulus.
    std::vector <double> ref_shear_moduli;
    for( auto& tri: triangles ){
        ref_shear_moduli.push_back(tri.ref_shear_modulus);
    }
    stuff.min_ref_shear_modulus = min_or_max_of_doubles(ref_shear_moduli, std::string("min"));

    stuff.char_flexural_rigidity = (1.0/6.0) * stuff.min_ref_shear_modulus * stuff.min_ref_thickness * stuff.min_ref_thickness * stuff.min_ref_thickness / (1.0-stuff.def_poisson_ratio);
    
    stuff.char_force_scale = stuff.min_ref_shear_modulus * stuff.approx_min_tri_size * stuff.min_ref_thickness;
    
    stuff.char_speed_scale = sqrt(stuff.min_ref_shear_modulus / stuff.ref_density);

    stuff.char_stretch_energy_density_scale = stuff.min_ref_shear_modulus * stuff.min_ref_thickness; // This is also a char membrane stress scale.
    // stuff.char_bend_energy_density_scale = stuff.min_ref_shear_modulus * pow(stuff.min_ref_thickness,3) / (stuff.SampleCharLength * stuff.SampleCharLength);
    
    std::vector <double> ref_areas;
    for( auto& tri: triangles ){
        ref_areas.push_back(tri.ref_area);
    }
    double tot_ref_area = kahan_sum(ref_areas);
    stuff.char_stretch_energy_scale = stuff.char_stretch_energy_density_scale * tot_ref_area;
    // stuff.char_bend_energy_scale = stuff.char_bend_energy_density_scale * totInitArea;
}


