/* 
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

Function to calculate a 3x6 matrix for each triangle that can be stored and
used repeatedly to obtain either the first or second spatial derivatives of the
patch for each triangle during the simulation. Forr each triangle t, the ids of 
the 3 out of 6 patch-defining nodes that are not vertices of t are also stored, 
as they will also be needed.
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

#include <cstddef>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <stdexcept>


#include "calc_deriv_mats.hpp"
#include "Tri_Class.hpp"
#include "Stuff_Class.hpp"
#include "Out_Stream_Class.hpp"


void calc_deriv_mats(
    Eigen::SparseMatrix<double, Eigen::RowMajor> &mat_for_continuum_quantities_from_dofs,
    Eigen::SparseMatrix<double, Eigen::RowMajor> &mat_for_continuum_quantities_from_dofs_transpose,
    const Eigen::Matrix<double,Eigen::Dynamic,3> &node_positions,
    const std::vector<Tri_Class> &triangles,
    const Stuff_Class &stuff,
    Out_Stream_Class &log){

    log << std::endl;

    // Using these funny things is just the recommended way of filling a sparse Eigen matrix.
    // Each triplet holds a row index, a column index, and the elements value.
    std::vector< Eigen::Triplet<double> >  mat_for_continuum_quantities_from_dofs_triplets;

    // dofs structure will be (node1_X, node2_X, ... , node1_Y, node2_Y, ... , node1_Z, node2_Z, ...)
    // continuum_quantities structure will be (tri1_Xx, tri1_Xy, tri1_Xxx, tri1_Xxy, tri1_Xyy, tri1_Yx, tri1_Yy, tri1_Yxx, tri1_Yxy, tri1_Yyy, tri1_Zx, tri1_Zy, tri1_Zxx, tri1_Zxy, tri1_Zyy, tri2_Xx, ...)
    
    
    for( const auto& tri: triangles ){

        // Index into continuum_quantities that is the start of triangle t's set of 15 continuum quantities
        int q = 15*tri.id;

        // Temp corresponds to finding the f = a + b(x-x_0) + c(y-y_0) that interpolates the values of some function f at tri's vertices.
        Eigen::Matrix<double,3,3> temp {
            {1.0, node_positions(tri.vertex_ids[0],0)-tri.ref_centroid(0), node_positions(tri.vertex_ids[0],1)-tri.ref_centroid(1)},
            {1.0, node_positions(tri.vertex_ids[1],0)-tri.ref_centroid(0), node_positions(tri.vertex_ids[1],1)-tri.ref_centroid(1)},
            {1.0, node_positions(tri.vertex_ids[2],0)-tri.ref_centroid(0), node_positions(tri.vertex_ids[2],1)-tri.ref_centroid(1)}
        };
        for(int v = 0; v < 3; ++v){

            // Xx
            mat_for_continuum_quantities_from_dofs_triplets.push_back(Eigen::Triplet<double>(q+0, tri.vertex_ids[v], temp.inverse()(1,v)));
            // Xy
            mat_for_continuum_quantities_from_dofs_triplets.push_back(Eigen::Triplet<double>(q+1, tri.vertex_ids[v], temp.inverse()(2,v)));

            // Yx
            mat_for_continuum_quantities_from_dofs_triplets.push_back(Eigen::Triplet<double>(q+5, tri.vertex_ids[v]+stuff.num_nodes, temp.inverse()(1,v)));
            // Yy
            mat_for_continuum_quantities_from_dofs_triplets.push_back(Eigen::Triplet<double>(q+6, tri.vertex_ids[v]+stuff.num_nodes, temp.inverse()(2,v)));

            // Zx
            mat_for_continuum_quantities_from_dofs_triplets.push_back(Eigen::Triplet<double>(q+10, tri.vertex_ids[v]+2*stuff.num_nodes, temp.inverse()(1,v)));
            // Zy
            mat_for_continuum_quantities_from_dofs_triplets.push_back(Eigen::Triplet<double>(q+11, tri.vertex_ids[v]+2*stuff.num_nodes, temp.inverse()(2,v)));
        }


        std::vector<int> patch_node_ids(tri.vertex_ids.begin(), tri.vertex_ids.end());
        patch_node_ids.insert(patch_node_ids.end(), tri.non_vertex_patch_nodes_ids.begin(), tri.non_vertex_patch_nodes_ids.end());
        for(int n = 0; n < 6; ++n){

            // Xxx
            mat_for_continuum_quantities_from_dofs_triplets.push_back(Eigen::Triplet<double>(q+2, patch_node_ids[n], tri.mat_for_patch_2nd_derivs(n,0)));
            // Xxy
            mat_for_continuum_quantities_from_dofs_triplets.push_back(Eigen::Triplet<double>(q+3, patch_node_ids[n], tri.mat_for_patch_2nd_derivs(n,1)));
            // Xyy
            mat_for_continuum_quantities_from_dofs_triplets.push_back(Eigen::Triplet<double>(q+4, patch_node_ids[n], tri.mat_for_patch_2nd_derivs(n,2)));

            // Yxx
            mat_for_continuum_quantities_from_dofs_triplets.push_back(Eigen::Triplet<double>(q+7, patch_node_ids[n]+stuff.num_nodes, tri.mat_for_patch_2nd_derivs(n,0)));
            // Yxy
            mat_for_continuum_quantities_from_dofs_triplets.push_back(Eigen::Triplet<double>(q+8, patch_node_ids[n]+stuff.num_nodes, tri.mat_for_patch_2nd_derivs(n,1)));
            // Yyy
            mat_for_continuum_quantities_from_dofs_triplets.push_back(Eigen::Triplet<double>(q+9, patch_node_ids[n]+stuff.num_nodes, tri.mat_for_patch_2nd_derivs(n,2)));

            // Zxx
            mat_for_continuum_quantities_from_dofs_triplets.push_back(Eigen::Triplet<double>(q+12, patch_node_ids[n]+2*stuff.num_nodes, tri.mat_for_patch_2nd_derivs(n,0)));
            // Zxy
            mat_for_continuum_quantities_from_dofs_triplets.push_back(Eigen::Triplet<double>(q+13, patch_node_ids[n]+2*stuff.num_nodes, tri.mat_for_patch_2nd_derivs(n,1)));
            // Zyy
            mat_for_continuum_quantities_from_dofs_triplets.push_back(Eigen::Triplet<double>(q+14, patch_node_ids[n]+2*stuff.num_nodes, tri.mat_for_patch_2nd_derivs(n,2)));
        }
    }

    mat_for_continuum_quantities_from_dofs.setFromTriplets(mat_for_continuum_quantities_from_dofs_triplets.begin(), mat_for_continuum_quantities_from_dofs_triplets.end());

    //mat_for_del_by_del_x.setFromTriplets(del_by_del_x_triplets.begin(), del_by_del_x_triplets.end());
    //mat_for_del_by_del_y.setFromTriplets(del_by_del_y_triplets.begin(), del_by_del_y_triplets.end());
    //mat_for_del_sq_by_del_x_sq.setFromTriplets(del_sq_by_del_x_sq_triplets.begin(), del_sq_by_del_x_sq_triplets.end());
    //mat_for_del_sq_by_del_x_del_y.setFromTriplets(del_sq_by_del_x_del_y_triplets.begin(), del_sq_by_del_x_del_y_triplets.end());
    //mat_for_del_sq_by_del_y_sq.setFromTriplets(del_sq_by_del_y_sq_triplets.begin(), del_sq_by_del_y_sq_triplets.end());

    mat_for_continuum_quantities_from_dofs.makeCompressed();
    //mat_for_del_by_del_x.makeCompressed();
    //mat_for_del_by_del_y.makeCompressed(); 
    //mat_for_del_sq_by_del_x_sq.makeCompressed(); 
    //mat_for_del_sq_by_del_x_del_y.makeCompressed(); 
    //mat_for_del_sq_by_del_y_sq.makeCompressed();

    mat_for_continuum_quantities_from_dofs_transpose = mat_for_continuum_quantities_from_dofs.transpose();

    //mat_for_del_by_del_x_transpose = mat_for_del_by_del_x.transpose();
    //mat_for_del_by_del_y_transpose = mat_for_del_by_del_y.transpose();
    //mat_for_del_sq_by_del_x_sq_transpose = mat_for_del_sq_by_del_x_sq.transpose();
    //mat_for_del_sq_by_del_x_del_y_transpose = mat_for_del_sq_by_del_x_del_y.transpose();
    //mat_for_del_sq_by_del_y_sq_transpose = mat_for_del_sq_by_del_y_sq.transpose();

    // dofs structure will be (node1_X, node2_X, ... , node1_Y, node2_Y, ... , node1_Z, node2_Z, ...)
    // continuum_quantities structure will be (tri1_Xx, tri1_Xy, tri1_Xxx, tri1_Xxy, tri1_Xyy, tri1_Yx, tri1_Yy, tri1_Yxx, tri1_Yxy, tri1_Yyy, tri1_Zx, tri1_Zy, tri1_Zxx, tri1_Zxy, tri1_Zyy, tri2_Xx, ...)
    //for( const auto& tri: triangles ){
    //    // matrix.block(i,j,p,q) means block of size (p,q) starting at (i,j).
    //    mat_for_continuum_quantities_from_dofs.block(tri.id, 0, 5, stuff.num_dofs) = mat_for_del_by_del_x.block(tri.id, 0, 1, stuff.num_dofs);
    //}
}
