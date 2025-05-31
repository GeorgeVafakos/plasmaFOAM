/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "voltageHandler.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
voltageHandler::voltageHandler(const fvMesh& mesh, const dictionary& dict, volScalarField& voltExt, const volScalarField& voltExtAmp) :
    mesh_(mesh),
    dict_(dict),
    voltExt_(voltExt),
    voltExtAmp_(voltExtAmp)
{ }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// AC voltage 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
voltageAC::voltageAC(const fvMesh& mesh, const dictionary& dict, volScalarField& voltExt, const volScalarField& voltExtAmp):
    voltageHandler(mesh, dict, voltExt, voltExtAmp)
{ }

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
volScalarField voltageAC::calcVoltage()
{
    scalar Vmax = dict_.getScalar("peakVoltage");
    scalar freq = dict_.getScalar("frequency");

    voltExt_ =  Vmax*voltExtAmp_*Foam::sin(2.0*M_PI*freq*mesh_.time().value());

    return voltExt_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Pulsed voltage 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
voltagePulsed::voltagePulsed(const fvMesh& mesh, const dictionary& dict, volScalarField& voltExt, const volScalarField& voltExtAmp):
    voltageHandler(mesh, dict, voltExt, voltExtAmp)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
volScalarField Foam::voltagePulsed::calcVoltage()
{
    // readParams();
    scalar Vmax = dict_.getScalar("peakVoltage");
    scalar freq = dict_.getScalar("frequency");
    scalar riseTime = dict_.getScalar("riseTime");
    scalar plateauWidth = dict_.getScalar("plateauWidth");
    scalar startPoint = dict_.lookupOrDefault<scalar>("startPoint", 0.0);
    scalar decayTime = dict_.lookupOrDefault<scalar>("decayTime", riseTime);
    scalar pulseWidth = 1.0/freq;

    try
    {
        if ( (startPoint+riseTime+plateauWidth+decayTime) >= pulseWidth )
        {
            throw std::runtime_error( "Error: The frequency value is lower than the sum of the pulse characteristic times!\n" );
            // Info<< "Error: The frequency value is lower than the sum of the pulse characteristic times!\n" << endl;
            // std::abort();
        }

        scalar numPeriods = floor(mesh_.time().value()/pulseWidth);
        scalar normalizedPulseTime = mesh_.time().value() - numPeriods*pulseWidth;

        voltExt_ =   Vmax*(voltExtAmp_/riseTime)*(normalizedPulseTime -startPoint) * ((normalizedPulseTime>=startPoint) && normalizedPulseTime<=startPoint+riseTime )  
                   + Vmax*(voltExtAmp_) * ((normalizedPulseTime>=startPoint+riseTime) && normalizedPulseTime<=startPoint+riseTime+plateauWidth )
                   + Vmax*(-voltExtAmp_/decayTime)*(normalizedPulseTime -startPoint-riseTime-plateauWidth-decayTime) * ((normalizedPulseTime>=startPoint+riseTime+plateauWidth) && normalizedPulseTime<=startPoint+riseTime+plateauWidth+decayTime )
                   + 0.0*voltExtAmp_ * ((normalizedPulseTime>=startPoint+riseTime+plateauWidth+decayTime) && normalizedPulseTime<=pulseWidth );

    }
    catch (const std::exception& e) 
    { 
        Pout<<e.what(); 
        std::abort();
    }

    return voltExt_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// DC voltage 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
voltageDC::voltageDC(const fvMesh& mesh, const dictionary& dict, volScalarField& voltExt, const volScalarField& voltExtAmp):
    voltageHandler(mesh, dict, voltExt, voltExtAmp)
{ }

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
volScalarField voltageDC::calcVoltage()
{
    scalar Vmax = dict_.getScalar("peakVoltage");

    voltExt_ =  Vmax*voltExtAmp_;

    return voltExt_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //





































// namespace Foam
// {

// voltageHandler::voltageHandler(const fvMesh& mesh, const dictionary& controlDict, const word& voltExtName, const word& voltExtAmpName) :
//     mesh_(mesh),
//     controlDict_(controlDict),
//     voltExtName_(voltExtName),
//     voltExtAmpName_(voltExtAmpName),
//     voltExt_(mesh.lookupObjectRef<volScalarField>(voltExtName)),
//     voltExtAmp_(mesh.lookupObjectRef<volScalarField>(voltExtAmpName))
// { }


// // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// // AC voltage 
// // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// voltageAC::voltageAC(const fvMesh& mesh, const dictionary& controlDict, const word& voltExtName, const word& voltExtAmpName):
//     voltageHandler(mesh, controlDict, voltExtName, voltExtAmpName)
// { }

// volScalarField voltageAC::calcVoltage()
// {
//     scalar Vmax = controlDict_.getScalar("peakVoltage");
//     scalar freq = controlDict_.getScalar("frequency");

//     voltExt_ =  Vmax*voltExtAmp_*Foam::sin(2.0*M_PI*freq*mesh_.time().value());

//     return voltExt_;
// }


// // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// // Pulsed voltage 
// // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// voltagePulsed::voltagePulsed(const fvMesh& mesh, const dictionary& controlDict, const word& voltExtName, const word& voltExtAmpName):
//     voltageHandler(mesh, controlDict, voltExtName, voltExtAmpName)
// {}

// volScalarField Foam::voltagePulsed::calcVoltage()
// {
//     // readParams();
//     scalar Vmax = controlDict_.getScalar("peakVoltage");
//     scalar freq = controlDict_.getScalar("frequency");
//     scalar riseTime = controlDict_.getScalar("riseTime");
//     scalar plateauWidth = controlDict_.getScalar("plateauWidth");
//     scalar startPoint = controlDict_.lookupOrDefault<scalar>("startPoint", 0.0);
//     scalar decayTime = controlDict_.lookupOrDefault<scalar>("decayTime", riseTime);
//     scalar pulseWidth = 1.0/freq;

//     try
//     {
//         if ( (startPoint+riseTime+plateauWidth+decayTime) >= pulseWidth )
//         {
//             throw std::runtime_error( "Error: The frequency value is lower than the sum of the pulse characteristic times!\n" );
//             // Info<< "Error: The frequency value is lower than the sum of the pulse characteristic times!\n" << endl;
//             // std::abort();
//         }

//         scalar numPeriods = floor(mesh_.time().value()/pulseWidth);
//         scalar normalizedPulseTime = mesh_.time().value() - numPeriods*pulseWidth;

//         voltExt_ =   Vmax*(voltExtAmp_/riseTime)*(normalizedPulseTime -startPoint) * ((normalizedPulseTime>=startPoint) && normalizedPulseTime<=startPoint+riseTime )  
//                     + Vmax*(voltExtAmp_) * ((normalizedPulseTime>=startPoint+riseTime) && normalizedPulseTime<=startPoint+riseTime+plateauWidth )
//                     + Vmax*(-voltExtAmp_/decayTime)*(normalizedPulseTime -startPoint-riseTime-plateauWidth-decayTime) * ((normalizedPulseTime>=startPoint+riseTime+plateauWidth) && normalizedPulseTime<=startPoint+riseTime+plateauWidth+decayTime )
//                     + 0.0*voltExtAmp_ * ((normalizedPulseTime>=startPoint+riseTime+plateauWidth+decayTime) && normalizedPulseTime<=pulseWidth );

//     }
//     catch (const std::exception& e) 
//     { 
//         Pout<<e.what(); 
//         std::abort();
//     }

//     return voltExt_;
// }

// // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// } // End namespace Foam

// // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
