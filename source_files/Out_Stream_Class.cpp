/* 
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

Non-templated member function definitions for Out_Stream_Class.
See Out_Stream_Class.hpp for details of the class.
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "Out_Stream_Class.hpp"


/* Set which file will be printed to (the log file, in the principal use
case). */
void Out_Stream_Class::set_output_filepath(const std::string &output_filepath) {

    filepath = output_filepath;
}


// Open or create log file corresponding to fileName member data.
void Out_Stream_Class::open(){

    // Open file in append mode.
    output_filestream.open(filepath, std::ofstream::app);

    if( !output_filestream ){
        throw std::runtime_error("Error: Problem creating or opening output_filestream.");
    }
}

// Close log file
void Out_Stream_Class::close(){
    output_filestream.close();
}

// Getter function to return a reference to the fileName.
const std::string& Out_Stream_Class::get_output_filepath() const {
    return filepath;
}

// Getter function to return a reference to the outputFileStream.
const std::ofstream& Out_Stream_Class::get_output_filestream() const {
    return output_filestream;
}


///////////////////////////////////////////////////////////////////////////

// For outputting manipulators such as std::endl.
/*
How this works: we first for convenience define a function pointer typedef so
that streamFunction is an alias for: "a pointer to a function that takes
a std::ostream& as input and returns a std::ostream&". std::endl is an
example of such a function. Then we overload the binary operator << to take
any such function, and use it as a SECOND argument. No first argument is
required because we're implementing this as a member function, so the first
argument is taken to be the calling object automatically!
Hence some_out_stream_class << std::endl takes some instance of Out_Stream_Class and feeds
that to its overloaded << operator, with the function (pointer I guess?)
std::endl as the second. What the << operator then does is apply the function
it took as an argument (std::endl here) and apply it in the way it usually
would be, first to std::cout, then to output_filestream. std::endl takes a stream
(output_filestream, say) and does it's newline-and-flush thing, and then returns
a reference to the stream so that you can then chain, i.e. put another << or
similar straight after. We want this behaviour too, which is why our
overloaded << operator has a return type and returns a reference to the
calling object (some_out_stream_class); so that chaining can occur.

See:
https://www.austincc.edu/comer/ds12s/binaryop.htm
Stackoverflow 4295432
Stackoverflow 3674200
Stackoverflow 14155364
*/

typedef std::ostream& (*stream_function) (std::ostream&);

Out_Stream_Class& Out_Stream_Class::operator << (const stream_function func){

    func(std::cout);
    open();
    func(output_filestream);
    close();
    return *this;
}

////////////////////////////////////////////////////////////////////////////////
