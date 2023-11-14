
//Turn Eigen bounds checking off.
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <Eigen/Dense>


#include "Node_Class.hpp"
#include "Tri_Class.hpp"
#include "Edge_Class.hpp"
#include "Stuff_Class.hpp"
#include "Out_Stream_Class.hpp"

void calc_mesh_info(
    std::vector<Node_Class> &nodes,
    const Eigen::Matrix<double,Eigen::Dynamic,3> &node_positions,
    std::vector<Tri_Class> &triangles,
    std::vector<Edge_Class> &edges,
    Stuff_Class &stuff,
    Out_Stream_Class &log);