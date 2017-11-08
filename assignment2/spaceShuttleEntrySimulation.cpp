/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

#include <Tudat/Astrodynamics/Aerodynamics/aerodynamicGuidance.h>

#include <assignment2/inputOutputDirectories.h>

namespace tudat
{


///                                                                 ////
/// Barebones class for aerodynamic guidance of the space shuttle   ////
///                                                                 ////
class SpaceShuttleAerodynamicGuidance: public aerodynamics::AerodynamicGuidance
{
public:

    //! Constructor
    SpaceShuttleAerodynamicGuidance( const simulation_setup::NamedBodyMap& bodyMap )
    {
        ////                                                                                ////
        ////    Retrieve relevant objects from environment and set as member variables      ////
        ////                                                                                ////
    }

    //! The aerodynamic angles are to be computed here
    void updateGuidance( const double time )
    {
        ////        currentAngleOfAttack_ = ...;        ////
        ////        currentAngleOfSideslip_ = ...;      ////
        ////        currentBankAngle_ = ...;            ////
    }

private:

};

}

//! Execute propagation of orbits of STS during entry.
int main( )
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    using namespace tudat::ephemerides;
    using namespace tudat::interpolators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::spice_interface;
    using namespace tudat::simulation_setup;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::propagators;
    using namespace tudat::aerodynamics;
    using namespace tudat::basic_mathematics;
    using namespace tudat::input_output;
    using namespace tudat;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ENVIRONMENT            //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels.
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );

    // Set simulation start epoch.
    const double simulationStartEpoch = 0.0;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ENVIRONMENT            //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define simulation body settings.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( { "Earth" } );
    bodySettings[ "Earth" ]->ephemerisSettings = boost::make_shared< simulation_setup::ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ), "SSB", "J2000" );
    bodySettings[ "Earth" ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    bodySettings[ "Earth" ]->atmosphereSettings = boost::make_shared< AtmosphereSettings >( nrlmsise00 );

    // Create Earth object
    simulation_setup::NamedBodyMap bodyMap = simulation_setup::createBodies( bodySettings );

    // Create vehicle objects.
    bodyMap[ "STS" ] = boost::make_shared< simulation_setup::Body >( );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Create vehicle  coefficients
    std::map< int, std::string > forceCoefficientFiles;
    forceCoefficientFiles[ 0 ] =
             tudat_applications::getStsInputPath( ) + "STS_CD.dat";
    forceCoefficientFiles[ 2 ] =
            tudat_applications::getStsInputPath( ) + "STS_CL.dat";

    boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            simulation_setup::readTabulatedAerodynamicCoefficientsFromFiles(
                forceCoefficientFiles, 2690.0 * 0.3048 * 0.3048,
                boost::assign::list_of( aerodynamics::angle_of_attack_dependent )(  aerodynamics::mach_number_dependent ),
                true, true );

    bodyMap[ "STS" ]->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "STS" ) );
    bodyMap[ "STS" ]->setConstantBodyMass( 165000 * 0.45 );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE ACCELERATIONS            ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define acceleration model settings.
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfSts;
    accelerationsOfSts[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfSts[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( aerodynamic ) );
    accelerationMap[  "STS" ] = accelerationsOfSts;

    bodiesToPropagate.push_back( "STS" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );


    ///                                                                         ///
    /// Create and set guidance model here                                      ///
    ///                                                                         ///

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             DEFINE INITIAL STATE                   ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set spherical elements for STS.
    Eigen::Vector6d stsSphericalEntryState;
    stsSphericalEntryState( SphericalOrbitalStateElementIndices::radiusIndex ) =
            spice_interface::getAverageRadius( "Earth" ) + 120.0E3;
    stsSphericalEntryState( SphericalOrbitalStateElementIndices::latitudeIndex ) = 0.3;
    stsSphericalEntryState( SphericalOrbitalStateElementIndices::longitudeIndex ) = 1.2;
    stsSphericalEntryState( SphericalOrbitalStateElementIndices::speedIndex ) = 7.45E3;
    stsSphericalEntryState( SphericalOrbitalStateElementIndices::flightPathIndex ) =
            -1.2 * mathematical_constants::PI / 180.0;
    stsSphericalEntryState( SphericalOrbitalStateElementIndices::headingAngleIndex ) = 0.6;

    // Convert sts state from spherical elements to Cartesian elements.
    Eigen::Vector6d systemInitialState = transformStateToGlobalFrame(
                convertSphericalOrbitalToCartesianState(
                                stsSphericalEntryState ),
                simulationStartEpoch, bodyMap.at( "Earth" )->getRotationalEphemeris( ) );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE INTEGRATION/PROPAGATION SETTINGS; PROPAGATE ORBIT; SAVE OUTPUT      ////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///                                                                                          ///
    ///    Create integration and propagation settings here, then propagate and save dynamics    ///
    ///                                                                                          ///

    std::string outputDirectory = tudat_applications::getStsOutputPath( );

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}

