/*
Function to write output data to output files:
a Legacy VTK PolyData file that can be read by ParaView
that contains mesh data plus things that you might 
want to "colour" the surface by, and an auxiliary 
file that holds some other things e.g. dial_factor.
Feel free to add things to either output, following
the same format.

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

#include <Eigen/Dense>

#include <cstddef>
#include <string>
#include <fstream>
#include <iomanip> //for setting output precision etc
#include <vector>
#include <stdexcept>

#include "write_output.hpp"
#include "Node_Class.hpp"
#include "Tri_Class.hpp"
#include "Stuff_Class.hpp"
#include "kahan_sum.hpp"

void write_output(
    const Eigen::Matrix<double,Eigen::Dynamic,1> &dofs,
    const std::vector<Tri_Class> &triangles,
    const int &stepcount,
    const double &dial_factor,
    const Eigen::Matrix<double,Eigen::Dynamic,2> &curvatures,
    const Eigen::Matrix<double,Eigen::Dynamic,2> &energy_densities,
    const Eigen::Matrix<double,Eigen::Dynamic,1> &angle_deficits,
    const Eigen::Matrix<double,Eigen::Dynamic,8> &stresses,
    const Eigen::Matrix<double,Eigen::Dynamic,1> &strains,
    [[maybe_unused]] const std::vector<double> &energies,
    const Stuff_Class &stuff){


    // First we write out the auxiliary data.
    std::ofstream aux_output;
    aux_output.open(stuff.aux_output_filepath, std::ofstream::app);
    if( stepcount == 0 ){
        aux_output << "stepcount  " <<  "dial_factor  " << "stuff.status  "
                    << "upper_slide_z  " << "upper_slide_force  " << "lower_slide_z  " << "lower_slide_force  " << std::endl; 
    }
    aux_output << stepcount << "  " <<  dial_factor << " " << stuff.status << "  " 
                << stuff.upper_slide_z_coord << " " << stuff.upper_slide_tot_vert_force << " " << stuff.lower_slide_z_coord << " " << stuff.lower_slide_tot_vert_force << " " << std::endl; 
    aux_output.close();

    // Now we write the vtk data.
    std::ofstream out_file(stuff.output_dir_path + "/stepcount_" + std::to_string(stepcount) + "_output.vtk");
    if( !out_file ){
        throw std::runtime_error("Error: Problem creating or opening vtk output file.");
    }

    out_file << std::scientific << std::setprecision(14)
    << "# vtk DataFile Version 4.2" << "\n"
    << "dial_factor = " << dial_factor << ", stepcount = " << stepcount;

    // Continue with rest of preamble and then output the relevant data to file.
    out_file << "\n" << "ASCII" << "\n"
    << "DATASET POLYDATA" << "\n"
    << "POINTS " << stuff.num_nodes << " double" << "\n";

    for(int n = 0; n < stuff.num_nodes; ++n){
        out_file << dofs(n) << " " << dofs(n+stuff.num_nodes) << " " << dofs(n+2*stuff.num_nodes) << "\n";
    }

    out_file << "POLYGONS " << stuff.num_tris << " " << 4*stuff.num_tris << "\n";

    for( const auto& tri: triangles ){
        out_file << "3 " << tri.vertex_ids[0] << " " << tri.vertex_ids[1] << " " << tri.vertex_ids[2] << "\n";
    }




    out_file << "CELL_DATA " << stuff.num_tris << "\n";

    out_file << "SCALARS gauss_curv double 1" << "\n"
    << "LOOKUP_TABLE default"<<"\n";

    for(int t = 0; t < stuff.num_tris; ++t){
        out_file << curvatures(t, 0) << "\n";
    }

    out_file << "SCALARS mean_curv double 1" << "\n"
    << "LOOKUP_TABLE default" << "\n";

    for(int t = 0; t < stuff.num_tris; ++t){
        out_file << curvatures(t, 1) << "\n";
    }

    out_file << "SCALARS dimless_stretch_energy_density double 1" << "\n"
    << "LOOKUP_TABLE default" << "\n";

    for(int t = 0; t < stuff.num_tris; ++t){
        out_file << energy_densities(t,0) / stuff.char_stretch_energy_density_scale << "\n";
    }

    out_file << "SCALARS dimless_bend_energy_density double 1" << "\n"
    << "LOOKUP_TABLE default" << "\n";

    for(int t = 0; t < stuff.num_tris; ++t){
        out_file << energy_densities(t,1) / stuff.char_stretch_energy_density_scale << "\n";
    }

    out_file << "SCALARS strain_measure double 1" << "\n"
    << "LOOKUP_TABLE default" << "\n";

    for(int t = 0; t < stuff.num_tris; ++t){
        out_file << strains(t) << "\n";
    }

    out_file << "FIELD cauchy_stress_info 4" << "\n";

    out_file << "dimless_cauchy_stress_eigval_1 1 " << stuff.num_tris << " double" << "\n";
    for(int t = 0; t < stuff.num_tris; ++t){
        out_file << stresses(t,0) / (triangles[t].ref_shear_modulus * triangles[t].ref_thickness) << "\n";
    }
    out_file << "cauchy_stress_eigvec_1 3 " << stuff.num_tris << " double" << "\n";
    for(int t = 0; t < stuff.num_tris; ++t){
        out_file << stresses(t,1) << " " << stresses(t,2) << " " << stresses(t,3) << "\n";
    }
    out_file << "dimless_cauchy_stress_eigval_2 1 " << stuff.num_tris << " double" << "\n";
    for(int t = 0; t < stuff.num_tris; ++t){
        out_file << stresses(t,4) / (triangles[t].ref_shear_modulus * triangles[t].ref_thickness) << "\n";
    }
    out_file << "cauchy_stress_eigvec_2 3 " << stuff.num_tris << " double" << "\n";
    for(int t = 0; t < stuff.num_tris; ++t){
        out_file << stresses(t,5) << " " << stresses(t,6) << " " << stresses(t,7) << "\n";
    }




    out_file << "POINT_DATA " << stuff.num_nodes << "\n";

    out_file << "SCALARS angleDeficit double 1" << "\n"
    << "LOOKUP_TABLE default" << "\n";

    for(int n = 0; n < stuff.num_nodes; ++n){
        out_file << angle_deficits(n) << "\n";
    }




    out_file.close();
}
