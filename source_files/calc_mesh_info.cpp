/*
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

Calculate a bunch of information about the mesh.

*/

//Turn Eigen bounds checking off.
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <stdexcept>
#include <cstddef> // For std::size_t


#include "calc_mesh_info.hpp"
#include "calc_triangles_incident_on_nodes.hpp"
#include "calc_node_neighbours.hpp"
#include "calc_triangle_adjacencies.hpp"
#include "calc_edges.hpp"
#include "Node_Class.hpp"
#include "Tri_Class.hpp"
#include "Edge_Class.hpp"
#include "Stuff_Class.hpp"
#include "Out_Stream_Class.hpp"
#include "kahan_sum.hpp"
#include "find_extremum_of_doubles.hpp"


void calc_mesh_info(
    std::vector<Node_Class> &nodes,
    const Eigen::Matrix<double,Eigen::Dynamic,3> &node_positions,
    std::vector<Tri_Class> &triangles,
    std::vector<Edge_Class> &edges,
    Stuff_Class &stuff,
    Out_Stream_Class &log){


    // Determine and store the ids of the triangles incident on each node.
    calc_triangles_incident_on_nodes( nodes, triangles);

    // Determine and store ids of each node's neighbours (other nodes that
    // share an edge with it).
    calc_node_neighbours( nodes, triangles );

    // Now determine, for each tri, which other tris it shares an edge with.
    calc_triangle_adjacencies( nodes, triangles);

    // Now calculate information about triangle edges.
    calc_edges( nodes, triangles, edges, stuff );

    stuff.num_dofs = 3 * stuff.num_nodes;
    stuff.num_continuum_quantities = 5 * 3 * stuff.num_tris; // 5 partial derivs for each of the 3 Cartesian components of 3D position, for each tri.


    // Set the edges debugging display functions to print to log.
    for( auto& edge: edges ){
        edge.edge_log_stream.set_output_filepath(log.get_output_filepath());
    }

    // Calculate and print the number of boundary and non-boundary edges.
    // Also id the nodes connected by boundary edges as boundary nodes.
    // Also calculate the total initial perimeter of the sample, to be used as a
    // characteristic sample length in estimating characteristic times, time steps
    // etc.
    stuff.num_boundary_edges = 0;
    std::vector<double> reference_boundary_edge_lengths(0);
    for( const auto& edge: edges ){
        if( edge.is_boundary == true ){

            stuff.num_boundary_edges += 1;

            nodes[edge.node_ids[0]].is_boundary = true;
            nodes[edge.node_ids[1]].is_boundary = true;

            // Store initial length of this boundary edge. 
            reference_boundary_edge_lengths.push_back( (node_positions.row(edge.node_ids[0]) - node_positions.row(edge.node_ids[1])).norm()  );
        }
    }
    // Calculate and store the number of boundary nodes.
    stuff.num_boundary_nodes = 0;
    for( const auto& node: nodes ){
        if( node.is_boundary == true ){
            stuff.num_boundary_nodes += 1;
        }
    }

    // Sum reference_boundary_edge_lengths to calculate perimeter, and 
    // take that to be the characteristic length scale of the sample's shape.
    double reference_perimeter = kahan_sum(reference_boundary_edge_lengths);
    stuff.sample_char_length = reference_perimeter;
    log << "Reference perimeter = " << reference_perimeter << std::endl;

    // Print some other things too.
    log << "Number of edges = " << stuff.num_edges << std::endl;
    log << "Number of boundary edges = " << stuff.num_boundary_edges << std::endl;
    log << "Number of non-boundary edges = " << stuff.num_edges - stuff.num_boundary_edges << std::endl;
    log << "Number of holes in mesh = " <<  1 + stuff.num_edges - stuff.num_nodes - stuff.num_tris  << std::endl; // From Euler's formula for a planar graph.
    log << "Number of boundary nodes = " << stuff.num_boundary_nodes << std::endl;

    // Further checks that things are ok:
    if( 3 * stuff.num_tris != 2 * stuff.num_edges - stuff.num_boundary_edges ){
        throw std::runtime_error(
            "Something has gone wrong in calculating triangle adjacencies and/or "
            "edges: the current edge and triangle counts violate a topological identity.");
    }
    if( 1 + stuff.num_edges - stuff.num_nodes - stuff.num_tris < 0 ){
        throw std::runtime_error(
            "Something is very wrong with the mesh, because the code thinks it has a "
            "negative number of holes! A first thing to check is that all nodes touch at least one tri.");
    }


    // Calculate and store each tri's reference area.
    std::vector<double> tri_ref_areas;
    for( auto& tri: triangles ){
        tri.calc_area(node_positions);
        tri.ref_area = tri.area;
        tri_ref_areas.push_back(tri.ref_area);
    }

    // Store the mesh's total reference area.
    stuff.total_ref_area = kahan_sum(tri_ref_areas);


    // Find the approximate smallest linear size of mesh element, based on
    // smallest altitude of each triangle, and also the worst aspect ratio,
    // which I define as smallest altitude / largest altitude.
    std::vector<double> smallest_altitudes;
    std::vector<double> aspect_ratios;
    for( const auto& tri: triangles ){

        Eigen::Vector3d side1 = node_positions.row(tri.vertex_ids[1]) - node_positions.row(tri.vertex_ids[0]);
        Eigen::Vector3d side2 = node_positions.row(tri.vertex_ids[2]) - node_positions.row(tri.vertex_ids[0]);
        Eigen::Vector3d side3 = node_positions.row(tri.vertex_ids[2]) - node_positions.row(tri.vertex_ids[1]);

        std::vector<double> altitudes {2.0 * tri.ref_area / side1.norm(),
                                       2.0 * tri.ref_area / side2.norm(),
                                       2.0 * tri.ref_area / side3.norm()};

        double smallest_alt = min_or_max_of_doubles(altitudes, std::string("min"));
        double largest_alt =  min_or_max_of_doubles(altitudes, std::string("max"));

        smallest_altitudes.push_back(smallest_alt);
        aspect_ratios.push_back(smallest_alt / largest_alt);
    }
    stuff.approx_min_tri_size = min_or_max_of_doubles(smallest_altitudes, std::string("min"));
    size_t smallest_tri = argmin_or_argmax_of_doubles(smallest_altitudes, std::string("min"));
    double worst_aspect_ratio = min_or_max_of_doubles(aspect_ratios, std::string("min"));
    size_t worst_aspect_tri = argmin_or_argmax_of_doubles(aspect_ratios, std::string("min"));

    log << "Approx smallest tri linear size = " << stuff.approx_min_tri_size << " (tri " << smallest_tri << ")" << std::endl;
    log << "Worst tri aspect ratio = " << worst_aspect_ratio << " (tri " << worst_aspect_tri << "). Around 0 is bad, 1 would be perfection." << std::endl;



    /* Calculate the cartesian coordinates in the initial flat state of each
    triangle's centroid, and store in a vector of 3-component vectors.*/
    for( auto& tri: triangles ){
        tri.ref_centroid = ( 
            node_positions.row(tri.vertex_ids[0]).transpose() + 
            node_positions.row(tri.vertex_ids[1]).transpose() + 
            node_positions.row(tri.vertex_ids[2]).transpose() )/3.0;
    }
    
}
