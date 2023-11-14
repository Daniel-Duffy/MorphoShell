/*
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

Function to calculate, for each triangle, a list of the triangles it shares
edges with, and store this list as member data.

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
#include <algorithm> // For std::sort etc.
#include <stdexcept>
#include <iostream>

#include "calc_triangle_adjacencies.hpp"
#include "Node_Class.hpp"
#include "Tri_Class.hpp"
#include "Edge_Class.hpp"
#include "Stuff_Class.hpp"

bool two_elements_in_common(std::array<int, 3> arr1, std::array<int, 3> arr2);

void calc_triangle_adjacencies(
    const std::vector<Node_Class> &nodes, 
    std::vector<Tri_Class> &triangles){

    // For each triangle (tri), ...
    for( auto& tri: triangles ){
        // ... we first make a list of all other
        // triangles that are incident on any of tri's vertices.
        std::vector<int> neighbour_tri_id_list(0);
        for( const auto& vert_id: tri.vertex_ids ){
            for( const auto& inc_tri_id: nodes.at(vert_id).incident_tri_ids ){
                if( inc_tri_id != tri.id ){
                    neighbour_tri_id_list.push_back(inc_tri_id);
                }
            }
        }

        // Then remove any duplicates in the list.
        std::sort(neighbour_tri_id_list.begin(), neighbour_tri_id_list.end());     
        auto last = std::unique(neighbour_tri_id_list.begin(), neighbour_tri_id_list.end());
        neighbour_tri_id_list.erase(last, neighbour_tri_id_list.end());

        // Then find which of the triangles in that list share *two* nodes with tri,
        // and therefore share an edge with tri. Then tri stores this info.
        // The two_elements_in_common function is defined down the page.
        tri.edge_sharing_tri_ids.resize(0);
        for( const auto& neighbour_tri_id: neighbour_tri_id_list ){
            if( two_elements_in_common(tri.vertex_ids, triangles[neighbour_tri_id].vertex_ids) ){
                tri.edge_sharing_tri_ids.push_back(neighbour_tri_id);
            }
        }        
    }


    // Finally do some checking of some simple things.
    // Note the first condition is not mathematically prohibited but would
    // indicate a very concerning mesh where some triangles share no edges with
    // other triangles.
    for( const auto& tri: triangles ){
        if(
           tri.edge_sharing_tri_ids.size() < 1
        || tri.edge_sharing_tri_ids.size() > 3
        ) {
            throw std::runtime_error(
                "Error: Problem with triangle (edge-sharing) adjacencies or edges setup. "
                "Either there is a bug, or some mesh pathology; investigate further!");
        }
    }
}




// Function to determine whether two three-element integer std::arrays have
// *exactly* two elements in common or not.
bool two_elements_in_common(std::array<int, 3> arr1, std::array<int, 3> arr2){

    std::sort(arr1.begin(), arr1.end());     
    std::sort(arr2.begin(), arr2.end()); 

    // Create vector to store common elements.
    std::vector<int> common;

    std::set_intersection(
        arr1.begin(), arr1.end(),
        arr2.begin(), arr2.end(),
        std::back_inserter(common));

    if( common.size() == 2 ){
        return true;
    }
    else{
        return false;
    }
}
