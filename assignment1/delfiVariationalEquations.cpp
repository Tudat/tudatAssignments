/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <Tudat/SimulationSetup/tudatEstimationHeader.h>

//! Get path for output directory.
static inline std::string getOutputPath(
        const std::string& extraDirectory = "" )
{
    // Declare file path string assigned to filePath.
    // __FILE__ only gives the absolute path of the header file!
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    std::string reducedPath = filePath_.substr( 0, filePath_.length( ) -
                                                std::string( "delfiVariationalEquations.cpp" ).length( ) );
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

//! Execute propagation of orbit of Delfi around the Earth.
int main()
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
    using namespace tudat::simulation_setup;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::gravitation;
    using namespace tudat::numerical_integrators;
    using namespace tudat::estimatable_parameters;
    using namespace tudat::ephemerides;


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation time settings.
    const double simulationStartEpoch = tudat::physical_constants::JULIAN_YEAR;
    const double simulationEndEpoch = simulationStartEpoch + 1.0 * tudat::physical_constants::JULIAN_DAY;

    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    // **** MODIFY FOR AE4867: define list of bodies required for simulation

    // Create body objects.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );

    // **** MODIFY FOR AE4867: update environment settings to comply with requirements

    NamedBodyMap bodyMap = createBodies( bodySettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create spacecraft object.
    bodyMap[ "Delfi" ] = boost::make_shared< simulation_setup::Body >( );

    bodyMap[ "Delfi" ]->setEphemeris( boost::make_shared< TabulatedCartesianEphemeris< > >(
                                            boost::shared_ptr< interpolators::OneDimensionalInterpolator
                                            < double, Eigen::Vector6d > >( ), "Earth", "ECLIPJ2000" ) );

    // **** MODIFY FOR AE4867: all required properties of body Delfi


    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // **** MODIFY FOR AE4867: define body that is propagated/central body

    // Define propagation settings.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfDelfi;
    basic_astrodynamics::AccelerationMap accelerationModelMap;

    // **** MODIFY FOR AE4867: define settings for accelerations acting on Delfi, and create acceleration models.



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set Keplerian elements for Delfi.

    // **** MODIFY FOR AE4867: set initial state according to your student number

    Eigen::Vector6d delfiInitialStateInKeplerianElements;
    delfiInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
    delfiInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
    delfiInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
//    delfiInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = 0.0;
//    delfiInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = 0.0;
//    delfiInitialStateInKeplerianElements( trueAnomalyIndex ) = 0.0;

    double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    const Eigen::Vector6d delfiInitialState = convertKeplerianToCartesianElements(
                delfiInitialStateInKeplerianElements, earthGravitationalParameter );


    // **** MODIFY FOR AE4867: create propagator and integrator settings.

    boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings;
    boost::shared_ptr< IntegratorSettings< > > integratorSettings;


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////    DEFINE PARAMETERS FOR WHICH SENSITIVITY IS TO BE COMPUTED   ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define list of parameters to estimate.
    std::vector< boost::shared_ptr< EstimatableParameterSettings > > parameterNames;

    // **** MODIFY FOR AE4867, question 3.2: define list of parameters for which state transition/sensitivity matrix are to be propagated.


    // Create parameters
    boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodyMap );

    // Print identifiers and indices of parameters to terminal.
    printEstimatableParameterEntries( parametersToEstimate );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT AND VARIATIONAL EQUATIONS         /////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // **** MODIFY FOR AE4867, question 3.1: propagate equations of motion
    // **** MODIFY FOR AE4867, question 3.2: propagate equations of motion and variational equations


    std::map< double, Eigen::MatrixXd > stateTransitionResult;
    std::map< double, Eigen::MatrixXd > sensitivityResult;
    std::map< double, Eigen::VectorXd > nominalIntegrationResult;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::string outputSubFolder = "DelfiVariationalEquations/";

    // Write perturbed satellite propagation history to file.
    input_output::writeDataMapToTextFile( nominalIntegrationResult,
                                          "singlePerturbedSatellitePropagationHistory.dat",
                                          getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( stateTransitionResult,
                                          "singlePerturbedSatelliteStateTransitionHistory.dat",
                                          getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( sensitivityResult,
                                          "singlePerturbedSatelliteSensitivityHistory.dat",
                                          getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    Eigen::Vector6d originalInitialState = delfiInitialState;

    // Iterate over all entries of initial state
    for( unsigned int entry = 0; entry < 6; entry++ )
    {
        // **** MODIFY FOR AE4867, question 3.3: propagate orbits needed to find limits of validity of linearization.
    }


    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}

