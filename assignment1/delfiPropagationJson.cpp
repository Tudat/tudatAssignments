/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <Tudat/JsonInterface/jsonInterface.h>

#include <boost/lexical_cast.hpp>

using namespace tudat::json_interface;
using namespace tudat::interpolators;
using namespace tudat;

//! Get path for output directory.
static inline std::string getOutputPath(
        const std::string& extraDirectory = "" )
{
    // Declare file path string assigned to filePath.
    // __FILE__ only gives the absolute path of the header file!
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    std::string reducedPath = filePath_.substr( 0, filePath_.length( ) -
                                                std::string( "delfiPropagationJson.cpp" ).length( ) );
    std::string outputPath = reducedPath + "SimulationOutput/";
    if( extraDirectory != "" )
    {
        outputPath += extraDirectory;
    }

    if( outputPath.at( outputPath.size( ) - 1 ) != '/' )
    {
        outputPath += "/";
    }

    return outputPath;
}

//! Execute propagation of orbits of Apollo during entry using the JSON Interface.
int main( )
{
    const std::string cppFilePath( __FILE__ );
    const std::string cppFolder = cppFilePath.substr( 0, cppFilePath.find_last_of("/\\") + 1 );

    JsonSimulationManager< > jsonSimulationManager( cppFolder + "delfiIntegratorAnalysis.json" );

    const std::string outputDirectory = getOutputPath( ) + "DelfiIntegrator/";

    // *********************** ITERATE OVER QUESTIONS ************************
    for( unsigned int question = 1; question <= 4; question++ )
    {
        jsonSimulationManager.resetJsonObject( jsonSimulationManager.getOriginalJsonObject( ) );

        unsigned int numberOfCases = 3;

        // *********************** ITERATE OVER ALL CASES FOR CURRENT QUESTION ************************
        for ( unsigned int i = 0; i < numberOfCases; ++i )
        {
            // Notify on propagation start
            std::cout << "Running propagation " << i + 1 << " of " << numberOfCases << " for question " <<question<< std::endl;

            // **** MODIFY FOR AE4867 question 2.2-2.4: update initial state if required

            if( question < 3 )
            {
                // **** MODIFY FOR AE4867, question 2.1 and 2.2: update settings for RK4 integrator and output files
                // for current question and case

                if( i == 0 )
                {                   
                    jsonSimulationManager[ "export" ][ 0 ][ "file" ] =
                            outputDirectory + "keplerOutputRk4_Case1_" + std::to_string( question ) + "_.dat";
                    jsonSimulationManager[ "export" ][ 1 ][ "file" ] =
                            outputDirectory + "cartesianOutputRk4_Case1_" + std::to_string( question ) + "_.dat";
                }
            }
            else
            {
                // **** MODIFY FOR AE4867, question 2.3 and 2.4: update settings for variable step-size integrator
                // and output files for current question and case

                jsonSimulationManager[ "integrator" ][ "type" ] = "rungeKuttaVariableStepSize";

            }

            // Create settings objects
            jsonSimulationManager.updateSettings( );

            // Propagate
            jsonSimulationManager.runPropagation( );

            // Export results
            jsonSimulationManager.exportResults( );

            // Silence unused key warnings after first propagation
            if ( i == 0 )
            {
                jsonSimulationManager[ "options" ][ "unusedKey" ] = tudat::json_interface::continueSilently;
            }

            // **** MODIFY FOR AE4867, question 2.4: interpolate results to common epochs

        }
    }

    return EXIT_SUCCESS;
}

