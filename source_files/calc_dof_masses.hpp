#include <vector>

#include "Node_Class.hpp"
#include "Tri_Class.hpp"
#include "Stuff_Class.hpp"
#include "Out_Stream_Class.hpp"


void calc_dof_masses(
    Eigen::Matrix<double,Eigen::Dynamic,1> &dof_masses,
    const std::vector<Tri_Class> &triangles,
    const Stuff_Class &stuff,
    [[maybe_unused]] Out_Stream_Class &log);