
#include <vector>
#include <cstddef>
#include <string>
#include <algorithm>
#include <stdexcept>

#include "find_extremum_of_doubles.hpp"

// Find minimum or maximum element in a std::vector of doubles.
double min_or_max_of_doubles(std::vector<double> &vec, const std::string &mode){

    std::vector<double>::iterator result;

    if( mode == "min" ){
        result = std::min_element(vec.begin(), vec.end());
    }
    else if( mode == "max" ){
        result = std::max_element(vec.begin(), vec.end());
    }
    else{
        throw std::runtime_error(
            "Error: provide a valid 2nd argument to min_or_max_of_doubles");
    }

    return *result;
}

// Find idx-of-minimum or idx-of-maximum in a std::vector of doubles.
size_t argmin_or_argmax_of_doubles(std::vector<double> &vec, const std::string &mode){

    std::vector<double>::iterator result;

    if( mode == "min" ){
        result = std::min_element(vec.begin(), vec.end());
    }
    else if( mode == "max" ){
        result = std::max_element(vec.begin(), vec.end());
    }
    else{
        throw std::runtime_error(
            "Error: provide a valid 2nd argument to min_or_max_of_doubles");
    }

    return std::distance(vec.begin(), result);
}