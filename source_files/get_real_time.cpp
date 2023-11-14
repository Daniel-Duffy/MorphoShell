/* 
/////////////////////////////////////////////////////
Copyright (C) 2020, Daniel Duffy, dld34@cam.ac.uk. All rights reserved.
Please cite Daniel Duffy and Dr John Biggins if you use any part of this 
code in work that you publish or distribute.

This file is part of Shellmorph.

Shellmorph is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Shellmorph is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Shellmorph.  If not, see <https://www.gnu.org/licenses/>.
/////////////////////////////////////////////////////

Function to calculate and return a std::string containing the current actual
date and time in some format. The format is just hardcoded in this function for
now.*/

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
