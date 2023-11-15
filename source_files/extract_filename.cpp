/*
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

Function to extract a filename from a full path, or just return the file name
if just a file name is given. The input is a std::string, and the function works
by taking just the piece of the string after the last / character, if any / is
present. */

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

#include <cstddef>
#include <string>

#include "extract_filename.hpp"

std::string extract_filename(const std::string &inputString){

    const size_t lastSlashLocation = inputString.find_last_of("/");

    //If no slash found, or if it is the last character just return the input
    //string
    if( lastSlashLocation == std::string::npos || lastSlashLocation == inputString.size()-1 ){
        return inputString;
    }
    //Otherwise return the part of the string after the last slash
    else{
        return inputString.substr(lastSlashLocation + 1);
    }
}
