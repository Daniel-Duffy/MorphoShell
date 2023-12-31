/*
This is the header file for the class that will be placed at each node,
containing the node's position, velocity etc.*/

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

#ifndef _NODE_CLASS_TAG_
#define _NODE_CLASS_TAG_

#include <Eigen/Dense>

#include "Out_Stream_Class.hpp"

class Node_Class{
public:

    /* Custom output stream allowing the debugging display function to print to
    a particular file in addition to std::cout.*/
    Out_Stream_Class node_log_stream;

    // id of this node, (also its index in the nodes container vector).
    int id;

    /* ids of the triangles with this node as a vertex ('incident triangles').
    Ordering is arbitrary.*/
    std::vector<int> incident_tri_ids;

    /* ids of this node's neighbours, i.e. those nodes connected to this node
    by triangle edges. Ordering is arbitrary.*/
    std::vector<int> neighbour_node_ids;

    /* Boolean representing whether the node is on a boundary of the sample
    (true) or not (false).*/
    bool is_boundary;

    // tag that will be read in from input file. Makes it easy to mark out certain
    // nodes for special treatment.
    bool tag;

    /* Bool representing whether the node's position is to be constrained
    in the impose_constraints function. This is just another kind of tag
    like the above.*/
    bool is_constrained;

    // Position vector (x, y and z coordinates).
    //Eigen::Vector3d pos;

    // Velocity vector.
    //Eigen::Vector3d vel;

    // Force vector.
    //Eigen::Vector3d force;

    // Mass assigned to node.
    double mass;

    /*Constructor, taking a single argument which is an output file name
    that gets the debugging display function to print to a particular file, as
    well as to std::out. This should usually be the log file (as for log).
    I ensure that default data values are recognisable values,
    for debugging.
    incident_tri_ids and neighbourNodeids are left with zero size at
    initialisation. */
    /*
    Node_Class(){
        id = -1234;
        seideIndicator = false;
        is_boundary = false;
        clamp_indicator = false;
        load_indicator = false;
        //pos.fill(-123456);
        vel.fill(654321);
        force.fill(56789);
        mass = -56789;
    }
    */

    // Declare other member functions.

    // Debugging function to display all member data.
    void display();

};
#endif
