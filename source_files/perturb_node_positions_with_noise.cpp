/*
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

Function to slightly perturb the node coordinates away from their current positions. 
The size of the disturbance is chosen to be small relative to the
approximate smallest element size, to ensure the right kind of scale.*/

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

#include <vector>
#include <Eigen/Dense>
#include <random>
#include <ctime>

#include "perturb_node_positions_with_noise.hpp"
#include "Node_Class.hpp"
#include "Stuff_Class.hpp"

void perturb_node_positions_with_noise(
    Eigen::Matrix<double,Eigen::Dynamic,3> &node_positions, 
    const Stuff_Class &stuff){

    // Set random number generator . A simple and common one is chosen here:
    // there is little point worrying about obtaining extremely 'good' random
    // numbers for such a simple task, in which the quality of randomness is
    // very unlikely to be important. The seed is chosen to be a fixed number,
    // rather than chosen based on time of day or anything like that, partly for 
    // simplicity, and partly because whether the time needs to be converted to 
    // an unsigned int or a long unsigned int seemed to be a bit delicate, and 
    // I think might make cross-platform portability more difficult.
    std::minstd_rand aSimpleEngine(123456);

    /// Set distribution to be symmetric about zero, and extend to a small
    // distance relative to the mesh triangle size. The 0.001 factor is hard coded but
    // could be changed if desired.
    std::uniform_real_distribution<double> distr(-stuff.approx_min_tri_size * 0.001, stuff.approx_min_tri_size * 0.001);

    for(int n = 0; n < stuff.num_nodes; ++n){
        for(int c = 0; c < 3; ++c){
            node_positions(n,c) += distr(aSimpleEngine);
        }
    }
}
