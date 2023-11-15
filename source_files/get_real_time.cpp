/* 
Function to calculate and return a std::string containing the current actual
date and time in some format. The format is just hardcoded in this function for
now.*/

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

#include <ctime>
#include <string>
#include <iomanip> 

#include "get_real_time.hpp"

std::string get_real_time(){

    std::time_t currRealTime = std::time(nullptr);
    std::tm realTimeStruct = *std::localtime(&currRealTime);
    std::stringstream realTimeStream;

    // Format chosen is day/month/year
    realTimeStream << std::put_time(&realTimeStruct,  "%d/%m/%y %H:%M");

    return realTimeStream.str();

}
