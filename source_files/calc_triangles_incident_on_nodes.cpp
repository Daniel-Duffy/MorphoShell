/*
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

Function to determine the ids of the triangles incident on each node (i.e.
having it as a vertex), and store these as node member data.*/

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

#include "calc_triangles_incident_on_nodes.hpp"
#include "Node_Class.hpp"
#include "Tri_Class.hpp"
#include "Stuff_Class.hpp"

void calc_triangles_incident_on_nodes(
    std::vector<Node_Class> &nodes, 
    const std::vector<Tri_Class> &triangles){

    // This first loop is just a precaution,
    // and almost certainly unnecessary.
    for( auto& node: nodes ){

        node.incident_tri_ids.resize(0);
    }

    for( const auto& tri: triangles ){

        // Add this triangle's id to each of its vertices in turn.
        for( const auto& vert: tri.vertex_ids ){
            nodes.at(vert).incident_tri_ids.push_back(tri.id);
        }
    }

}
