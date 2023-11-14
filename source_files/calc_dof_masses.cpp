/* 
/////////////////////////////////////////////////////
Copyright (C) 2020, Daniel Duffy, dld34@cam.ac.uk. All rights reserved.
Please cite Daniel Duffy and Dr John Biggins if you use any part of this 
code in work that you publish or distribute.

This file is part of Shellmorph.

Shellmorph is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Shellmorph is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Shellmorph.  If not, see <https://www.gnu.org/licenses/>.
/////////////////////////////////////////////////////

Function to finish setting up initial data for the flat LCE sheet: setting
node velocities to zero, and storing the initial in-plane sides' components for
the triangles (using the x-y plane basis), which will then not change. The
initial areas are also calculated and stored. Also, calculate node masses by
having each triangle contribute 1/3 of its initial mass to each of its vertcies.
*/

//Turn Eigen bounds checking off for speed.
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <vector>
#include <Eigen/Dense>
#include <cmath>

#include "calc_dof_masses.hpp"
#include "Node_Class.hpp"
#include "Tri_Class.hpp"
#include "Stuff_Class.hpp"
#include "Out_Stream_Class.hpp"

void calc_dof_masses(
    Eigen::Matrix<double,Eigen::Dynamic,1> &dof_masses,
    const std::vector<Tri_Class> &triangles,
    const Stuff_Class &stuff,
    [[maybe_unused]] Out_Stream_Class &log){


    for(int d = 0; d < stuff.num_dofs; ++d){
        //Set all nodes masses to zero before calculating them next
        dof_masses(d) = 0.0;
    }

    for( const auto& tri: triangles ){

        /* Add 1/3 of the mass of each triangle to each of its vertices.*/
        for(int v = 0; v < 3; ++v){
            dof_masses[tri.vertex_ids[v]] += stuff.ref_density * tri.ref_area * tri.ref_thickness / 3.0;
            dof_masses[tri.vertex_ids[v]+stuff.num_nodes] += stuff.ref_density * tri.ref_area * tri.ref_thickness / 3.0;
            dof_masses[tri.vertex_ids[v]+2*stuff.num_nodes] += stuff.ref_density * tri.ref_area * tri.ref_thickness / 3.0;
        }
    }
}




















