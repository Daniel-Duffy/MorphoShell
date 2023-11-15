/* 
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

This is based on the code given in an answer to Stackoverflow 14155364.
See Out_Stream_Class.cpp for function details.

Class to be used like std::cout << "text here" << std::endl etc, replacing
std::cout with an instance of Out_Stream_Class that you've constructed:
someInstance.

The Out_Stream_Class constructor used to construct someInstance should be
passed the path and name of a text file. Then using :
someInstance <<
instead of:
std::cout <<
will write whatever follows BOTH to std::cout AND to the chosen file. This is
useful if you want to write a log file but also want everything in the log
file to get written to std::cout e.g. for easy viewing in a terminal.
*/

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

#ifndef _OUT_STREAM_CLASS_TAG_
#define _OUT_STREAM_CLASS_TAG_

#include <iostream>
#include <fstream>
#include <string>

class Out_Stream_Class
{

public:

    // Constructor (just default).
    Out_Stream_Class() {};

    /* Set path+name of file which will be printed to (e.g. a log file, in the
    principal use case). */
    void set_output_filepath(const std::string&);

    // Open or create output file corresponding to fileName member data.
    void open();

    // Close output file.
    void close();

    // Getter function to return a reference to the fileName.
    const std::string& get_output_filepath() const;

    // Getter function to return a reference to the outputFileStream.
    const std::ofstream& get_output_filestream() const;

    /* For regular output of variables etc.
    This is defined here since the function is templated so defining it just in
    the .cpp file leads to problems with compilation/linking!*/
    template<typename T> Out_Stream_Class& operator << (const T& something){

    std::cout << something;
    open();
    output_filestream << something;
    close();

    return *this;
    }

    ///////////////////////////////////////////////////////////////////////////
    // For outputting manipulators such as std::endl.
    typedef std::ostream& (*stream_function) (std::ostream&);
    Out_Stream_Class& operator << (const stream_function);
    //////////////////////////////////////////////////////////////////////////

    private:
        // File path+name the file stream will correspond to.
        std::string filepath;

        // Actual output file stream
        std::ofstream output_filestream;

};

#endif
