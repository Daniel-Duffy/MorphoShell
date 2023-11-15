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
#include <omp.h>


#include "calc_angle_deficits.hpp"
#include "Node_Class.hpp"
#include "Tri_Class.hpp"
#include "Stuff_Class.hpp"
#include "kahan_sum.hpp"

void calc_angle_deficits(
    Eigen::Matrix<double,Eigen::Dynamic,1> &angle_deficits,
    const Eigen::Matrix<double,Eigen::Dynamic,1> &dofs,
    const std::vector<Node_Class> &nodes,
    const std::vector<Tri_Class> &triangles,
    [[maybe_unused]] const Stuff_Class &stuff){

    const double pi = 3.14159265358979323846;

    #pragma omp parallel for simd
    for( int n = 0; n < stuff.num_nodes; ++n ){

        if( nodes[n].is_boundary ){
            // I think this is the natural choice 
            // given that there is an exact discrete
            // version of Gauss-Bonnet for a mesh of
            // triangles.
            angle_deficits(n) = pi;
            
        }
        else{
            angle_deficits(n) = 2.0 * pi;
        }
    }

    // Do NOT multi-thread loops like this where different
    // threads would try to write to the same location in 
    // memory at the same time!
    for( const auto& tri: triangles ){

        Eigen::Vector3d vert_0 {dofs(tri.vertex_ids[0]), dofs(tri.vertex_ids[0]+stuff.num_nodes), dofs(tri.vertex_ids[0]+2*stuff.num_nodes)};
        Eigen::Vector3d vert_1 {dofs(tri.vertex_ids[1]), dofs(tri.vertex_ids[1]+stuff.num_nodes), dofs(tri.vertex_ids[1]+2*stuff.num_nodes)};
        Eigen::Vector3d vert_2 {dofs(tri.vertex_ids[2]), dofs(tri.vertex_ids[2]+stuff.num_nodes), dofs(tri.vertex_ids[2]+2*stuff.num_nodes)};

        Eigen::Vector3d side_0 = vert_1 - vert_0;
        Eigen::Vector3d side_1 = vert_2 - vert_1;
        Eigen::Vector3d side_2 = -side_1 - side_0;

        double side_length_0 = side_0.norm();
        double side_length_1 = side_1.norm();
        double side_length_2 = side_2.norm();

        angle_deficits[tri.vertex_ids[0]] += -acos( -(side_2.dot(side_0)) / (side_length_2 * side_length_0) );

        angle_deficits[tri.vertex_ids[1]] += -acos( -(side_0.dot(side_1)) / (side_length_0 * side_length_1) );

        angle_deficits[tri.vertex_ids[2]] += -acos( -(side_1.dot(side_2)) / (side_length_1 * side_length_2) );
    }
}
