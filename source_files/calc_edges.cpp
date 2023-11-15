/*
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

Function to calculate all manner of things to do with
mesh edges in the reference state.
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
#include <vector>
#include <cstddef> // For std::size_t
#include <limits> // For INT_MAX
#include <algorithm> // For std::find
#include <stdexcept>

#include "calc_edges.hpp"
#include "Node_Class.hpp"
#include "Tri_Class.hpp"
#include "Edge_Class.hpp"
#include "Stuff_Class.hpp"

void calc_edges(
    const std::vector<Node_Class> &nodes, 
    std::vector<Tri_Class> &triangles, 
    std::vector<Edge_Class> &edges, 
    Stuff_Class &stuff){

    // We don't yet know how many edges there are, but it will be a lot most of
    // the time, so we reserve suitable space using an upper bound on the
    // number of edges. We'll get rid of any excess space later.
    edges.reserve(3*stuff.num_tris);

    // This is a vector in which each element corresponds to the node
    // at the same position in the nodes container. If the element equals
    // 1, that means that the corresponding node has already added its 
    // edges to the edges container, while a value of 0 means it hasn't.
    // I use 0s and 1s instead of bools because std::vector<bool> is
    // very strange (not a normal std::vector) and is to be avoided.
    std::vector<int> already_added_edges_to_container(stuff.num_nodes, 0);

    // For each node, ...
    for( const auto& node: nodes ){

        // we look at each of its neighbour (edge-sharing) nodes
        // If a neighbour has already accounted for the shared edge, 
        // we move on, but if not we add the shared edge to the edges container.
        for( const auto& neighbour_id: node.neighbour_node_ids ){
            if( already_added_edges_to_container.at(neighbour_id) == 1 ){
                continue;
            }
            else{
                edges.emplace_back();
                edges.back().node_ids[0] = node.id;
                edges.back().node_ids[1]  = neighbour_id;

                if(edges.size() > INT_MAX){
                    throw std::overflow_error("The unsigned quantity edges.size() overflowed upon conversion to int (i.e. was larger than INT_MAX)");
                }
                else{
                    edges.back().id = static_cast<int>(edges.size()) - 1;
                }
            }
        }

        // Now that this node has added its edges to the edges
        // container, we record that so those edges won't be
        // double-counted later in the loop.
        already_added_edges_to_container.at(node.id) = 1;
    }


    // Now shrink capacity of std::vector edges to match its actual number of
    // elements (edges) (i.e. just allow discarding of extra unused storage).
    edges.shrink_to_fit();

    // Store total number of edges.
    if(edges.size() > INT_MAX)
    {
        throw std::overflow_error("The unsigned quantity edges.size() overflowed upon conversion to int (i.e. was larger than INT_MAX)");
    }
    else{
        stuff.num_edges = static_cast<int>(edges.size());
    }


    // Now we work out which triangles each edge is an edge of. 
    // To do so, we use the fact that those triangles are exactly
    // the (either 1 or 2) triangles that are incident on both of
    // an edge's nodes.
    for( auto& edge: edges ){

        std::vector<int> node0_inc_tris = nodes[edge.node_ids[0]].incident_tri_ids;
        std::vector<int> node1_inc_tris = nodes[edge.node_ids[1]].incident_tri_ids;
        std::sort(node0_inc_tris.begin(), node0_inc_tris.end());     
        std::sort(node1_inc_tris.begin(), node1_inc_tris.end()); 
        // Create vector to store elements common to both lists of incident tri ids.
        std::vector<int> common;
        std::set_intersection(
            node0_inc_tris.begin(), node0_inc_tris.end(),
            node1_inc_tris.begin(), node1_inc_tris.end(),
            std::back_inserter(common));

        if( common.size() == 1){
            edge.is_boundary = true;
            edge.adjacent_tri_ids.resize(1);
            edge.adjacent_tri_ids[0] = common[0];
        }
        else if( common.size() == 2 ){
            edge.is_boundary = false;
            edge.adjacent_tri_ids.resize(2);
            edge.adjacent_tri_ids[0] = common[0];
            edge.adjacent_tri_ids[1] = common[1];
        }
        else{
            throw std::runtime_error(
        "Something has gone very wrong in calc_edges, because the nodes at each end of an edge "
        "should have either one or two incident triangles in common, but the code thinks it's "
        "found an exception. Investigate!");
        }

    }


    // Now we make sure that triangles know which edges in the 
    // edges container are their edges.
    std::vector<int> cursor_idxs_into_tri_edge_ids(stuff.num_tris, 0);
    for( const auto& edge: edges ){
        for( auto& adj_tri_id: edge.adjacent_tri_ids ){
            triangles.at(adj_tri_id).edge_ids[cursor_idxs_into_tri_edge_ids[adj_tri_id]] = edge.id;
            cursor_idxs_into_tri_edge_ids[adj_tri_id] += 1;
        }
    }
}
