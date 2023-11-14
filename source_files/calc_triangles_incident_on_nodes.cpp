/*
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

Function to determine the ids of the triangles incident on each node (i.e.
having it as a vertex), and store these as node member data.*/

//Turn Eigen bounds checking off for speed.
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
