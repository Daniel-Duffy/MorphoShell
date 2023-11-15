/* 
This is the header file for the class that will correspond to each edge,
containing the ids of the nodes on the edge and at the other corners of the
two triangles sharing the edge, and the ids of those triangles.*/

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

#ifndef _EDGE_CLASS_TAG_
#define _EDGE_CLASS_TAG_

#include <Eigen/Dense>

#include "Out_Stream_Class.hpp"

class Edge_Class{
public:

    /* Custom output stream allowing the debugging display function to print to
    a particular file in addition to std::cout.*/
    Out_Stream_Class edge_log_stream;

    // id so this edge 'knows' which it is
    int id;

    /* ids (and indexes in the nodes container vector) of the nodes that the
    edge is defined to start and end at.*/
    std::array<int, 2> node_ids;

    /* ids (and indices in the triangles container vector) of the (either 1 or
    2) triangles that this edge is an edge of. We term these triangles 'adjacent'
    to the edge.*/
    std::vector<int> adjacent_tri_ids;

    /* Boolean representing whether the edge is on the boundary of the sample
    (true) or not (false).*/
    bool is_boundary;

    /*Constructor, taking a single argument which is an output file name
    that gets the debugging display function to print to a particular file, as
    well as to std::out. This should usually be the log file (as for log).
    I ensure that default data values are recognisable values,
    for debugging. */
    Edge_Class(){
        id = -54321;
        //node_ids.fill(-1234);
        is_boundary = false;
    }

    // Declare other member functions.

    // Debugging function to display all member data.
    void display();

};
#endif
