/*
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

Function to extract a filename from a full path, or just return the file name
if just a file name is given. The input is a std::string, and the function works
by taking just the piece of the string after the last / character, if any / is
present. */

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
