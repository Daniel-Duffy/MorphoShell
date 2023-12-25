/* 
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

Function to calculate a 3x6 matrix for each triangle that can be stored and
used repeatedly to obtain either the first or second spatial derivatives of the
patch for each triangle during the simulation. For each triangle t, the ids of 
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
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <string> // For std::to_string

#include "set_up_patch_fitting.hpp"
#include "Node_Class.hpp"
#include "Tri_Class.hpp"
#include "Stuff_Class.hpp"
#include "Out_Stream_Class.hpp"


void set_up_patch_fitting(
    const std::vector<Node_Class> &nodes,
    const Eigen::Matrix<double,Eigen::Dynamic,3> &node_positions,
    std::vector<Tri_Class> &triangles,
    const Stuff_Class &stuff,
    Out_Stream_Class &log){

    int num_non_boundary_tris_that_tried_multiple_patch_choices = 0;

    for( auto& tri: triangles ){

        /* We now take a simple approach of assigning the
        non_vertex_patch_nodes_ids to nodes which are not vertices of the given
        triangle, but are the closest nodes to the given triangle's initial
        (flat state) centroid, (chosen from the surrounding triangles' vertices),
        that give a reasonable answer for the matForPatchFirstDerivs. Some
        choices of nodes lead to the inversion of a near-singular matrix en
        route to matForPatchFirstDerivs, giving poor results, hence the
        'reasonable' caveat. The common problematic patch node configurations
        are two parallel rows of three nodes each, or anything where more than 3 
        nodes lie along a single straight line. If all the patch node choice possibilites
        explored here are exhausted without success even for moderate condition number thresholds,
        the search could be extended.*/

        // Loop over vertices of current triangle.
        std::vector<double> dists_to_centroid;
        std::vector<int> possible_patch_node_ids;
        for(int v = 0; v < 3; ++v){

            // Loop over the triangles incident on each of these vertices.
            for(size_t j = 0; j < nodes[tri.vertex_ids[v]].incident_tri_ids.size(); ++j){

                // Loop over the vertices of these incident triangles.
                for(int w = 0; w < 3; ++w){

                    int this_node_id = triangles[nodes[tri.vertex_ids[v]].incident_tri_ids[j]].vertex_ids[w];
                    bool is_node_already_accounted_for = false;

                    /* Now check whether this node is a vertex of the
                    triangle under consideration.*/
                    for(int u = 0; u < 3; ++u){
                        if( this_node_id == tri.vertex_ids[u] ){
                            is_node_already_accounted_for = true;
                        }
                    }

                    /* Check whether current node has already been added to
                    list of potential patch nodes.*/
                    for(size_t u = 0; u < possible_patch_node_ids.size(); ++u){
                        if( this_node_id == possible_patch_node_ids[u] ){
                            is_node_already_accounted_for = true;
                        }
                    }

                    /* If the current node is not a vertex, and is not already in
                    the list of potentials, find its distance to the current
                    triangle's initial centroid position, and store the
                    distance to this node along with the others in a list of
                    potentials. Store also this node's identity id similarly.
                    */
                    if( is_node_already_accounted_for == false ){
                        dists_to_centroid.push_back( (node_positions.row(this_node_id).transpose() - tri.ref_centroid).norm() );
                        possible_patch_node_ids.push_back( this_node_id );
                    }
                }
            }
        }


        /* Now get a vector listing the *indices* (in the dists_to_centroid
        vector) of distances to the current triangle centroid, sorted from
        smallest to largest distance. This can then be used to extract the
        corresponding possible_patch_node_ids in order of distance. */
        std::vector<size_t> idx_in_dist_list(dists_to_centroid.size());
        for(size_t q = 0; q < idx_in_dist_list.size(); ++q){
               idx_in_dist_list[q] = q;
        }
        std::sort(std::begin(idx_in_dist_list), std::end(idx_in_dist_list),[&dists_to_centroid](const size_t &idx1, const size_t &idx2)->bool{ return dists_to_centroid[idx1] < dists_to_centroid[idx2]; });
        // idx_in_dist_list is now sorted by distance to centroid.

        // Put closest node in non_vertex_patch_nodes_ids.
        tri.non_vertex_patch_nodes_ids[0] = possible_patch_node_ids[idx_in_dist_list[0]];

        bool success = false;
        
        // Now loop over some candidates for non_vertex_patch_nodes_ids[1]. So far I explore
        // just the three next-closest nodes. But that search could be extended.
        for(std::size_t p = 1; p < std::min(static_cast<std::size_t>(3), idx_in_dist_list.size()); ++p){
            tri.non_vertex_patch_nodes_ids[1] = possible_patch_node_ids.at(idx_in_dist_list.at(p));
        

            /* Now loop over possibilities for non_vertex_patch_nodes_ids[2], finding
            the first choice that gives reasonable results where the matrix to be
            inverted is not close to singular.*/

            for(size_t q = 2; q < possible_patch_node_ids.size(); ++q){

                tri.non_vertex_patch_nodes_ids[2] = possible_patch_node_ids[idx_in_dist_list[q]];


                /* Calculate an approx `linear size' for the patch, taken to be the
                RMS distances of the 6 patch nodes to the central triangle's
                centroid. This measure is chosen to minimise the influence of any
                one node, and, to be roughly unaffected by the geometric considerations
                that lead to a badly conditioned problem. This leads to the chosen
                threshold test having roughly the same level of `strictness' for each
                triangle.*/
                double this_patch_size = sqrt((
                    (node_positions.row(tri.vertex_ids[0]).transpose() - tri.ref_centroid).squaredNorm()+
                    (node_positions.row(tri.vertex_ids[1]).transpose() - tri.ref_centroid).squaredNorm()+
                    (node_positions.row(tri.vertex_ids[2]).transpose() - tri.ref_centroid).squaredNorm()+
                    (node_positions.row(tri.non_vertex_patch_nodes_ids[0]).transpose() - tri.ref_centroid).squaredNorm()+
                    (node_positions.row(tri.non_vertex_patch_nodes_ids[1]).transpose() - tri.ref_centroid).squaredNorm()+
                    (node_positions.row(tri.non_vertex_patch_nodes_ids[2]).transpose() - tri.ref_centroid).squaredNorm()
                    ) / 6.0);

                /* Calculate a matrix for each patch, that when multiplied
                onto the current patch node coords, gives the coefficients specifying the
                quadratic surface that goes through all 6 nodes in the patch. We will then
                only retain and store the bottom three rows, which give the second
                derivatives of position over the patch (wrt 2D parametrisation coords),
                which is all that's needed for our 2nd F.F. approximation. If the full
                quadratic surface was needed in future, the full matrix could easily be
                retained and used instead.
                If the 3D deformed state coordinates as a function of the 2D parametrisation
                coordinates are X(x,y), Y(x,y), Z(x,y), we have the approximations:
                X(x,y) = a_1 + b_1*(x-x_0) + c_1*(y-y_0) + d_1*(x-x_0)^2 + e_1*(y-y_0)^2 + f_1*(x-x_0)(y-y_0)
                Y(x,y) = a_2 + b_2*(x-x_0) + ..... and equivalently for Z(x,y) etc,
                where we approximate the surface as quadratic in the vicinity of a triangle
                centroid at x_0, y_0. We do this for each triangle, and use the known
                (x,y) and (X,Y,Z) coordinates of each of the patch nodes as constraints on
                the coefficients a_1, b_1,....,f_3. Finding these coefficients then
                simply requires inverting a 6x6 matrix, as below.*/

                // Loop over patch nodes for this triangle.
                Eigen::Matrix<double,6,6> temp_patch_node_data_matrix;
                for(int n = 0; n < 6; ++n){

                    /* Patch nodes will always stick to this order: the three vertex
                    nodes first, then the others.*/
                    Eigen::Vector3d patch_node_pos;
                    if( n < 3 ){
                        patch_node_pos = node_positions.row(tri.vertex_ids[n]).transpose();
                    }
                    else{
                        patch_node_pos = node_positions.row(tri.non_vertex_patch_nodes_ids[n-3]).transpose();
                    }

                    /* NB the following assumes that the initial state is a flat sheet,
                    the (x,y) coords on which then are the 2D parametrisation
                    coordinates, which then distort with the deforming sheet. If the
                    initial state were not flat, we would need to pick some less
                    obvious parametrisation to construct this matrix, and to define the
                    deformation gradient with respect to. The factors of 0.5 ensure
                    means that down the line, if we apply a part of the inverse of
                    this matrix, what we get are second derivatives of the patch,
                    rather than coefficients of the quadratic terms in the patch
                    expansion.*/
                    temp_patch_node_data_matrix(0,n) = 1;
                    temp_patch_node_data_matrix(1,n) = (patch_node_pos(0) - tri.ref_centroid(0));
                    temp_patch_node_data_matrix(2,n) = (patch_node_pos(1) - tri.ref_centroid(1));
                    temp_patch_node_data_matrix(3,n) = 0.5 * temp_patch_node_data_matrix(1,n) * temp_patch_node_data_matrix(1,n);
                    temp_patch_node_data_matrix(4,n) = temp_patch_node_data_matrix(1,n) * temp_patch_node_data_matrix(2,n);
                    temp_patch_node_data_matrix(5,n) = 0.5 * temp_patch_node_data_matrix(2,n) * temp_patch_node_data_matrix(2,n);
                }

                /* For this triangle:
                (a_1, b_1, c_1, d_1, e_1, f_1) * temp_patch_node_data_matrix= (X1, X2, X3, X4, X5, X6)
                where X1, X2,..,X6 will be the X coordinates of the 6 patch nodes in any
                3D deformed state. Similarly for the 'Y's and 'Z's with (a_2, b_2,...)
                and (a_3, b_3,...) respectively. Thus inverting temp_patch_node_data_matrix
                allows us to find the 'a's, 'b's,..., 'f's from any deformed state
                3D coords of the 6 patch nodes. We check first that
                temp_patch_node_data_matrix is actually invertible.
                */
                bool mat_is_invertible = true;
                Eigen::Matrix<double,6,3> temp_mat_for_patch_2nd_derivs;
                double temp_abs_cond_num;
                double temp_abs_cond_num_divided_by_thresh_val;
                Eigen::FullPivLU< Eigen::Matrix<double,6,6> > temp_LU(temp_patch_node_data_matrix);
                if( !temp_LU.isInvertible() ){
                    mat_is_invertible = false;
                }
                else{
                    temp_mat_for_patch_2nd_derivs = temp_patch_node_data_matrix.inverse().block<6,3>(0,3);


                    /* Compute absolute condition number i.e. the maximum singular
                    value of the matrix which will be used to find the patch
                    derivatives, to check that it will give reasonable results. If a
                    certain direction is not well 'sampled' by the patch for example,
                    the condition number will be high, and the results of using it
                    would be poor (and can lead to overly strict time step
                    restrictions), so an alternative choice of patch nodes should be
                    considered instead.*/
                    Eigen::JacobiSVD<Eigen::Matrix<double,6,3> > temp_SVD(temp_mat_for_patch_2nd_derivs);
                    temp_abs_cond_num = temp_SVD.singularValues()(0); //This has dimensions 1/Length^2.
                    temp_abs_cond_num_divided_by_thresh_val = temp_abs_cond_num * this_patch_size * this_patch_size / stuff.patch_mat_dimless_conditioning_thresh;
                }

                /* Only proceed if this matrix will give reasonable results, i.e. has
                a low enough condition number. The threshold can be tuned in the
                settings file. It will be mesh dependent, with poorer meshes
                giving higher condition numbers, and thus requiring a higher
                threshold to actually find any patch that is considered good enough.
                The search may be extended in future versions of the code to help
                with this, if it becomes necessary.*/
                if( (!(temp_abs_cond_num_divided_by_thresh_val < 1.0)) || mat_is_invertible == false ){

                    /* The first time the above condition occurs for a triangle,
                    add one to a counter if the triangle in question is *not* on the
                    boundary. This is because if lots of interior triangles are
                    going through multiple patch possibilities before satisfying the
                    determinant threshold conditions,
                    stuff.PatchMatrixDimensionlessConditioningThreshold is likely
                    too high, or something may be strange about the mesh. */
                    if( q==2 && !tri.is_boundary ){
                        num_non_boundary_tris_that_tried_multiple_patch_choices += 1;
                    }

                    //Move onto the next possible non-vertex patch node since this
                    //one has produced an unsuitable matrix.
                    continue;
                }
                else{
                    /* Store the matrix that gives the first or second derivatives
                    of the patch, depending on which secFF approx is being used.*/

                    tri.mat_for_patch_2nd_derivs = temp_mat_for_patch_2nd_derivs;

                    //Move on, as search has now been successfully
                    //completed.
                    success = true;
                    break;
                }
            }
            if( success ){
                break;
            }
        }
        //Throw error to main if whole search has been exhausted unsuccessfully for this tri.
        if( !success ) {
            throw std::runtime_error("At least one search for patch nodes was exhausted without "
            "success (triangle " + std::to_string(tri.id) + "); all possible patch matrices in the search had a  "
            "condition number above the acceptance threshold. Try increasing this threshold. If "
            "that does not solve the issue, or causes other issues, either the patch node search will "
            "need to be extended, or you will need to use a nicer mesh. Aborting.");
        }

    }

     /*Print also the number of non-boundary triangles that had to search through
     multiple possible patch options to find one satisfying the determinant
     condition. It seems unlikely that the determinant will be a small for
     triangles in the interior of a reasonable mesh, so if this number is large,
     that suggests something suspicious. One explanation might be that
     stuff.PatchMatrixDimensionlessConditioningThreshold has been set to too
     low a value*/
     log << "\nNumber of non-boundary triangles that had to search through \n" <<
     "multiple possible patch options to find one \nsatisfying the condition number " <<
     "criterion was " << num_non_boundary_tris_that_tried_multiple_patch_choices <<
     ", \nwhich should not be a large proportion of the mesh's triangles." << std::endl;
}
