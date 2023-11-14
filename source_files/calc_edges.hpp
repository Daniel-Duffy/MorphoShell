#include <vector>

#include "Node_Class.hpp"
#include "Tri_Class.hpp"
#include "Edge_Class.hpp"
#include "Stuff_Class.hpp"

void calc_edges(
    const std::vector<Node_Class> &nodes, 
    std::vector<Tri_Class> &triangles, 
    std::vector<Edge_Class> &edges, 
    Stuff_Class &stuff);