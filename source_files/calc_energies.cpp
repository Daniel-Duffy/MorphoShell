/*
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

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

#include <vector>
#include <Eigen/Dense>

#include "calc_energies.hpp"
#include "Node_Class.hpp"
#include "Tri_Class.hpp"
#include "Stuff_Class.hpp"
#include "kahan_sum.hpp"

void calc_energies(
    std::vector<double> &energies,
    const Eigen::Matrix<double,Eigen::Dynamic,2> &energy_densities,
    const Eigen::Matrix<double,Eigen::Dynamic,1> &velocities,
    const Eigen::Matrix<double,Eigen::Dynamic,1> &dof_masses,
    const std::vector< Tri_Class > &triangles,
    const Stuff_Class &stuff){

    std::vector<double> temp(stuff.num_tris);
    #pragma omp parallel for simd
    for(int t = 0; t < stuff.num_tris; ++t){
        temp[t] = triangles[t].ref_area * energy_densities(t,0);
    }
    energies[0] = kahan_sum(temp); // Total stretch energy.

    #pragma omp parallel for simd
    for(int t = 0; t < stuff.num_tris; ++t){
        temp[t] = triangles[t].ref_area * energy_densities(t,1);
    }
    energies[1] = kahan_sum(temp); // Total bend energy.

    temp.resize(stuff.num_nodes);
    #pragma omp parallel for simd
    for(int n = 0; n < stuff.num_nodes; ++n){
        Eigen::Vector3d vel_vec {velocities(n), velocities(n+stuff.num_nodes), velocities(n+2*stuff.num_nodes)};
        temp[n] = 0.5 * dof_masses(n) * vel_vec.squaredNorm();
    }
    energies[2] = kahan_sum(temp); // Total kinetic energy.

}
