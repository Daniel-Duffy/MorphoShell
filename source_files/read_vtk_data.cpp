/*
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

Function that takes a data filename and reads node coordinates and the
node ids of each triangle into prepared std:vectors of the node and triangle
classes. The number of nodes and triangles is also put into the settings struct.
The programmed metric tensors and second fundamental forms are then also read
into the appropriate data structures.
The file is assumed to correspond to a text file containing VTK Legacy
PolyData (triangles specifically). See VTK documentation for explanation of
format. NODE LABELS START AT ZERO. In fact the assumed format is even more
stringent than just vtk (vtk reader would be required otherwise); it must be
exactly as in the example data files accompanying this code.

If a filename is supplied as a third command line argument, a dial-in factor,
counter into the programmed tensor sequence, and a set of node positions are read
in from that file. This is to provide the option of beginning evolution from an
ansatz, which could for example be where a previous unfinished simulation left
off. The format for such ansatz files must therefore be exactly as is used for
the OUTPUT vtk files, except that only the initial preamble (minus time and
stepcount), and the node position data are required.*/

//For testing (spaces can be fiddly with ignore()!) can use:
/*
std::string teststring;
ref_data_file >> teststring;
std::cout << teststring << std::endl;
*/

//Turn Eigen bounds checking off.
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <iostream>
#include <cstddef>
#include <fstream>
#include <string>
#include <vector>
#include <limits>
#include <Eigen/Dense>
#include <stdexcept>

#include "read_vtk_data.hpp"
#include "Node_Class.hpp"
#include "Tri_Class.hpp"
#include "Stuff_Class.hpp"
#include "Out_Stream_Class.hpp"

void read_vtk_data(
    std::vector<Node_Class> &nodes,
    Eigen::Matrix<double,Eigen::Dynamic,3> &node_positions,
    std::vector< Eigen::Vector3d > &ansatz_node_positions,
    std::vector<Tri_Class> &triangles,
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &abar_info,
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &bbar_info,
    Stuff_Class &stuff,
    double &dial_factor_to_start_from,
    Out_Stream_Class &log){


    // Temporary string variable used to help check that 
    // the data files have the correct format. 
    std::string temp;

    log << "Beginning reading of reference data file." << std::endl;


    std::ifstream ref_data_file(stuff.ref_data_filepath);
    if( !ref_data_file ){
        throw std::runtime_error("Error: Problem opening reference data file.");
    }

    // Ignore first 4 lines.
    for(int i=0; i<4; ++i){
        ref_data_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    // Ignore "POINTS".
    ref_data_file.ignore(std::numeric_limits<std::streamsize>::max(), ' ');

    // Get number of nodes and resize nodes containers accordingly.
    ref_data_file >> stuff.num_nodes;
    if( !(stuff.num_nodes > 0) ){
        throw std::runtime_error("Error: Problem with (or before) line giving number of nodes in reference data file.");
    }
    else{
        nodes.resize(stuff.num_nodes);
        node_positions.resize(stuff.num_nodes,3);
    }
    // Set the nodes debugging display functions to print to the log file.
    for( auto& node: nodes  ){
        node.node_log_stream.set_output_filepath(log.get_output_filepath());
    }


    // Ignore rest of line.
    ref_data_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // Put node id and coordinates in nodes container.
    for( int node_id = 0; node_id < stuff.num_nodes; ++node_id )
    {
        nodes.at(node_id).id = node_id;
        ref_data_file >> node_positions(node_id, 0);
        ref_data_file >> node_positions(node_id, 1);
        ref_data_file >> node_positions(node_id, 2);
        ref_data_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Ignore delimeter at end of line.
    }


    // Check we've arrived at the expected line of the file: 'POLYGONS'
    ref_data_file >> temp;
    if( temp != "POLYGONS" ){
        throw std::runtime_error("Error: Problem with reference data file, at or before 'POLYGONS' line. E.g. may have provided or stated wrong number of nodes.");
    }
    temp.clear();


    /* Get number of triangles and resize triangles container. */
    ref_data_file >> stuff.num_tris;
    if( !(stuff.num_tris > 0) ){
        throw std::runtime_error("Error: Problem with line giving number of triangles in reference data file.");
    }
    else{
        triangles.resize(stuff.num_tris);
    }
    // Set the triangles debugging display functions to print to the log file.
    for( auto& tri: triangles  ){
        tri.tri_log_stream.set_output_filepath(log.get_output_filepath());
    }

    // Ignore rest of line.
    ref_data_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // Put triangles' self and vertex ids in triangles container.
    for( int tri_id = 0; tri_id < stuff.num_tris; ++tri_id )
    {
        triangles.at(tri_id).id = tri_id;
        // Ignore first number of each row, which just says this polygon has
        // 3 vertices.
        ref_data_file.ignore(1);
        ref_data_file >> triangles.at(tri_id).vertex_ids[0];
        ref_data_file >> triangles.at(tri_id).vertex_ids[1];
        ref_data_file >> triangles.at(tri_id).vertex_ids[2];
        ref_data_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    // Check that the triangulation is using correct node id convention.
    bool found_zero = false;
    for( const auto& tri: triangles ){
        for( int n = 0; n < 3; ++n ){
            if( tri.vertex_ids[n] == 0 ){
                found_zero = true;
            }
        }
    } 
    if( not found_zero ){
        throw std::runtime_error(
            "Error: MorphoShell requires the triangulation you provide to use node ids starting from 0. "
            "It appears that the triangulation you provided did not do this, since it didn't contain a 0 anywhere.");
    }


    // Check we've arrived at the expected line of the file: 'CELL_DATA'...
    ref_data_file >> temp;
    if( temp != "CELL_DATA" ){
        throw std::runtime_error(
            "Error: Problem with reference data file, at or before 'CELL_DATA' line. "
            "E.g. you may have provided or stated wrong number of triangles.");
    }
    temp.clear();
    // Ignore rest of line.
    ref_data_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // Ignore line starting FIELD...
    ref_data_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // Check we've arrived at the expected line of the file: 'abar_info'...
    ref_data_file >> temp;
    if( temp != "abar_info" ){
        throw std::runtime_error("Error: Problem with reference data file, at or before an abar_info line.");
    }
    temp.clear();
    // Read how many pieces of abar_info each tri has, and resize
    // the storage container accordingly.
    int num_pieces_of_abar_info_per_tri;
    ref_data_file >> num_pieces_of_abar_info_per_tri;
    abar_info.resize(stuff.num_tris, num_pieces_of_abar_info_per_tri);
    // Ignore rest of line.
    ref_data_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    // Read abar_info.
    for( int tri_id = 0; tri_id < stuff.num_tris; ++tri_id ){
        for( int p = 0; p < num_pieces_of_abar_info_per_tri; ++p ){

        ref_data_file >> abar_info(tri_id, p);
        }
        ref_data_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    // Check we've arrived at the expected line of the file: bbar_info...
    ref_data_file >> temp;
    if( temp != "bbar_info" ){
        throw std::runtime_error(
            "Error: Problem with reference data file, at or before a bbar_info line. "
            "E.g. you may have provided wrong number of abar_info lines.");
    }
    temp.clear();
    // Read how many pieces of bbar_info each tri has, and resize
    // the storage container accordingly.
    int num_pieces_of_bbar_info_per_tri;
    ref_data_file >> num_pieces_of_bbar_info_per_tri;
    bbar_info.resize(stuff.num_tris, num_pieces_of_bbar_info_per_tri);
    // Ignore rest of line.
    ref_data_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    // Read bbar_info.
    for( int tri_id = 0; tri_id < stuff.num_tris; ++tri_id ){
        for( int p = 0; p < num_pieces_of_bbar_info_per_tri; ++p ){

        ref_data_file >> bbar_info(tri_id, p);
        }
        ref_data_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }


    // Check we've arrived at the expected line of the file: ref_thicknesses...
    ref_data_file >> temp;
    if( temp != "ref_thicknesses" ){
        throw std::runtime_error(
            "Error: Problem with reference data file, at or before a ref_thicknesses line. "
            "E.g. you may have provided wrong number of bbar_info lines.");
    }
    temp.clear();
    // Ignore rest of line.
    ref_data_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    // Put triangles' ref_thickness values in their container.
    if( stuff.ref_thickness_if_uniform > 0 ){
        log << "Using uniform ref_thickness = " << stuff.ref_thickness_if_uniform << std::endl;
    }
    else{
        log << "Using non-uniform ref_thickness from vtk file. " << std::endl;
    }
    for( int tri_id = 0; tri_id < stuff.num_tris; ++tri_id )
    {
        ref_data_file >> triangles.at(tri_id).ref_thickness;
        ref_data_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        // Override with uniform value from settings file if a positive uniform
        // value was provided.
        if( stuff.ref_thickness_if_uniform > 0 ){
            triangles.at(tri_id).ref_thickness = stuff.ref_thickness_if_uniform;
        }
    }
    
    // Check we've arrived at the expected line of the file: ref_shear_moduli...
    ref_data_file >> temp;
    if( temp != "ref_shear_moduli" ){
        throw std::runtime_error(
            "Error: Problem with reference data file, at or before a ref_shear_moduli line. "
            "E.g. you may have provided wrong number of ref_thicknesses lines.");
    }
    temp.clear();
    // Ignore rest of line.
    ref_data_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    // Put triangles' ref_shear_moduli values in their container.
    if( stuff.ref_shear_modulus_if_uniform > 0 ){
        log << "Using uniform ref_shear_modulus = " << stuff.ref_shear_modulus_if_uniform << std::endl;
    }
    else{
        log << "Using non-uniform ref_shear_modulus from vtk file. " << std::endl;
    }
    for( int tri_id = 0; tri_id < stuff.num_tris; ++tri_id )
    {
        ref_data_file >> triangles.at(tri_id).ref_shear_modulus;
        ref_data_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        // Override with uniform value from settings file if a positive uniform
        // value was provided.
        if( stuff.ref_shear_modulus_if_uniform > 0 ){
            triangles.at(tri_id).ref_shear_modulus = stuff.ref_shear_modulus_if_uniform;
        }
    }

    // Check we've arrived at the expected line of the file: tri_tags...
    ref_data_file >> temp;
    if( temp != "tri_tags" ){
        throw std::runtime_error("Error: Problem with reference data file, at or before tri_tags.");
    }
    temp.clear();
    // Ignore rest of line.
    ref_data_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    // Read in tri tags.
    // std::ios_base::boolalpha ensures that '0's and '1's  will be interpreted as bools.
    ref_data_file.unsetf(std::ios_base::boolalpha);
    for( int tri_id = 0; tri_id < stuff.num_tris; ++tri_id )
    {
        ref_data_file >> triangles.at(tri_id).tag;
        ref_data_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    // Check we've arrived at the expected line of the file: 'POINT_DATA'...
    ref_data_file >> temp;
    if( temp != "POINT_DATA" ){
        throw std::runtime_error(
            "Error: Problem with reference data file, at or before 'POINT_DATA'. "
            "E.g. you may have provided wrong number of ref_thickness values.");
    }
    temp.clear();
    // Ignore rest of line.
    ref_data_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    // Ignore next line.
    ref_data_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // Check we've arrived at the expected line of the file: constraint_indicators...
    ref_data_file >> temp;
    if( temp != "constraint_indicators" ){
        throw std::runtime_error("Error: Problem with reference data file, at or before constraint_indicators.");
    }
    temp.clear();
    // Ignore rest of line.
    ref_data_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    // Read in node constraint indicators.    
    // std::ios_base::boolalpha ensures that '0's and '1's  will be interpreted as bools.
    ref_data_file.unsetf(std::ios_base::boolalpha);
    for( int node_id = 0; node_id < stuff.num_nodes; ++node_id )
    {
        ref_data_file >> nodes.at(node_id).is_constrained;
        ref_data_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    // Check we've arrived at the expected line of the file: node_tags...
    ref_data_file >> temp;
    if( temp != "node_tags" ){
        throw std::runtime_error("Error: Problem with reference data file, at or before node_tags.");
    }
    temp.clear();
    // Ignore rest of line.
    ref_data_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    // Read in node tags.
    // std::ios_base::boolalpha ensures that '0's and '1's  will be interpreted as bools.
    ref_data_file.unsetf(std::ios_base::boolalpha);
    for( int node_id = 0; node_id < stuff.num_nodes; ++node_id )
    {
        ref_data_file >> nodes.at(node_id).tag;
        ref_data_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    /* Close the data file, as we should be at the end now. We check we
    actually did reach the end first, and for other errors that can be caught
    automatically by a std::ifstream. */
    if( ref_data_file.bad() || ref_data_file.fail() ){
        throw std::runtime_error(
            "Error: Unkown problem with reference data file. One of .bad(), .fail() was true - investigate! "
            "E.g. you may not have supplied enough node indicator data in the file, so the file may have reached its end prematurely.");
    }
    ref_data_file >> temp; // Attempts to read whatever's left, which should be nothing!
    if( !ref_data_file.eof() ){
        throw std::runtime_error("Error: Did not reach end of reference data file as expected. Check that it has exactly the correct format.");
    }
    ref_data_file.close();

    log << "Completed reading of reference data file. Note; it is still possible that something went wrong in the process.\n" << std::endl;

    //////////////////////////////////////////////////////////////////////
    // Now read the ansatz file if there is one. 
    // Only the preamble and node position data is read;
    // whatever comes after that is ignored.

    if( stuff.ansatz_file_was_given ){

        log << "Beginning reading of ansatz data file." << std::endl;

        std::ifstream ansatz_DataFile(stuff.ansatz_filepath);
        if(!ansatz_DataFile){
            throw std::runtime_error("Error: Problem opening ansatz data file.");
        }


        // Ignore first line.
        ansatz_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');


        /* Get the dial factor that evolution of the ansatz should start from.*/
        ansatz_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
        ansatz_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
        ansatz_DataFile >> dial_factor_to_start_from;
        // Ignore rest of line.
        ansatz_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        // Check a few things that could have gone wrong.
        if( !(0 <= dial_factor_to_start_from) || !(dial_factor_to_start_from<= 1) ){
            throw std::runtime_error("Error: Problem reading ansatz data file. Remember it must have exactly the correct format.");
        }

        // Ignore two lines of preamble.
        ansatz_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        ansatz_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        // Check number of nodes stated in ansatz file matches stuff.num_nodes.
        ansatz_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
        int ansatz_stated_num_nodes;
        ansatz_DataFile >> ansatz_stated_num_nodes;
        if( !(ansatz_stated_num_nodes == stuff.num_nodes) ){
            throw std::runtime_error(
                "Error: Problem with (or before) line giving number of nodes in ansatz data file. "
                "E.g. you may have supplied ansatz and ref data files with inconsistent numbers of nodes.");
        }
        // Ignore rest of line.
        ansatz_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');


        // Put node ansatz coordinates into suitable container.
        ansatz_node_positions.resize(stuff.num_nodes);
        for( int node_id = 0; node_id < stuff.num_nodes; ++node_id )
        {
            ansatz_DataFile >> ansatz_node_positions.at(node_id)(0);
            ansatz_DataFile >> ansatz_node_positions.at(node_id)(1);
            ansatz_DataFile >> ansatz_node_positions.at(node_id)(2);
            ansatz_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }


        /* Close the data file, ignoring any remaining data. We also check
        for other errors that can be caught automatically by a
        std::ifstream. */
        if( ansatz_DataFile.bad() || ansatz_DataFile.fail() ){
            throw std::runtime_error(
                "Error: Unkown problem with ansatz data file. One of .bad(), .fail() was true - investigate! "
                " E.g. you may not have supplied enough ansatz node positions, so the file reached its end prematurely.");
        }
        ansatz_DataFile.close();

        log << "Completed reading of ansatz data file. Note, it is still possible that something went wrong in the process.\n" << std::endl;
    }
}
