/*
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

Function to calculate, for each triangle, a list of the triangles it shares
edges with, and store this list as member data, as well as a list of edge ids
for the triangle. For a given triangle, the n elements of the edge-sharing
triangles list correspond to the first n elements of the edge ids list, with
the remaining edge ids corresponding to boundary edges.

The edges data structure is also set up, where each edge stores its two end
nodes, its adjacent triangles etc.

There may well be a more efficient approach to that taken here, with less vector resizing, 
but this only occurs once and is therefore unlikely to be a bottleneck worth worrying about.
*/

//Turn Eigen bounds checking off for speed.
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <Eigen/Dense>
#include <vector>
#include <cstddef> // For std::size_t
#include <algorithm> // For std::find
#include <stdexcept>

#include "calc_node_neighbours.hpp"
#include "Node_Class.hpp"
#include "Tri_Class.hpp"
#include "Stuff_Class.hpp"

void calc_node_neighbours(
    std::vector<Node_Class> &nodes, 
    const std::vector<Tri_Class> &triangles){

    // For each node, ...
    for( auto& node: nodes ){
        
        // ... we make a list of all the
        // other nodes that are vertices of 
        // triangles incident on the node.
        // This set of nodes is sometimes called
        // a 1-ring.
        node.neighbour_node_ids.resize(0); // Precautionary.
        for( const auto& inc_tri_id: node.incident_tri_ids ){
            for( const auto& vert: triangles.at(inc_tri_id).vertex_ids ){
                if( vert != node.id ){
                    node.neighbour_node_ids.push_back(vert);
                }
            }
        }

        // Then we remove any duplicates in the list.
        std::sort(node.neighbour_node_ids.begin(), node.neighbour_node_ids.end());     
        auto last = std::unique(node.neighbour_node_ids.begin(), node.neighbour_node_ids.end());
        node.neighbour_node_ids.erase(last, node.neighbour_node_ids.end());
    }
}
