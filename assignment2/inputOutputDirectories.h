#ifndef TUDAT_ASSIGNMENT_INPUTOUTPUTDIRECTORIES_H
#define TUDAT_ASSIGNMENT_INPUTOUTPUTDIRECTORIES_H

#include <string>

namespace tudat_applications
{

//! Get path for output directory.
static inline std::string getStsOutputPath( )
{
    // Declare file path string assigned to filePath.
    // __FILE__ only gives the absolute path of the header file!
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    return ( filePath_.substr( 0, filePath_.length( ) -
                             std::string( "inputOutputDirectories.h" ).length( ) ) ) + "StsResults/";
}

//! Get path for output directory.
static inline std::string getStsInputPath( )
{
    // Declare file path string assigned to filePath.
    // __FILE__ only gives the absolute path of the header file!
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    return ( filePath_.substr( 0, filePath_.length( ) -
                             std::string( "inputOutputDirectories.h" ).length( ) ) ) + "StsInput/";
}

}

#endif // INPUTOUTPUTDIRECTORIES_H
