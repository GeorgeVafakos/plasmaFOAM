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

#include "transportCoeffsHandler.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(transportCoeffsHandler, 0);
defineRunTimeSelectionTable(transportCoeffsHandler, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
transportCoeffsHandler::transportCoeffsHandler(const fvMesh& mesh, const dictionary& dict, volScalarField& transpCoeff) :
    mesh_(mesh),
    dict_(dict),
    transpCoeff_(transpCoeff),
    E_(mesh.lookupObjectRef<volVectorField>("E")),
    N_(mesh.lookupObjectRef<volScalarField>("NGas")),
    P_(mesh.lookupObjectRef<volScalarField>("P")),
    EtoN_("EtoN", 1.0e21*mag(E_/N_)), 
    EtoP_("EtoP", mag(E_/P_)), 
    dimEtoN_("dimEtoN", EtoN_.dimensions(), scalar(1)),
    dimEtoP_("dimEtoP", EtoP_.dimensions(), scalar(1)),
    dimMu_("dimMu", dimensionSet(-1, 0, 2, 0, 0, 1, 0), scalar(1)),
    dimD_("dimD", dimViscosity, scalar(1)),
    dimN_("dimN", N_.dimensions(), scalar(1)),
    dimP_("dimP", P_.dimensions(), scalar(1)),
    ones_
    (
        IOobject
        (
            "ones",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("temp",dimless,scalar(1))
    )
{ }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Zero coefficients
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
transpCoeffsNone::transpCoeffsNone(const fvMesh& mesh, const dictionary& dict, volScalarField& transpCoeff):
    transportCoeffsHandler(mesh, dict, transpCoeff)
{ }

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void transpCoeffsNone::calcMobilityCoeffs()
{
    transpCoeff_ = dimensionedScalar("temp", dimensionSet(-1, 0, 2, 0, 0, 1, 0), 0.0);
}

void transpCoeffsNone::calcDiffusionCoeffs()
{
    transpCoeff_ = dimensionedScalar("temp", dimViscosity, 0.0);
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Constant coefficients
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
transpCoeffsConstValue::transpCoeffsConstValue(const fvMesh& mesh, const dictionary& dict, volScalarField& transpCoeff):
    transportCoeffsHandler(mesh, dict, transpCoeff)
{ }

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void transpCoeffsConstValue::calcMobilityCoeffs()
{
    scalar coeffValue(dict_.get<scalar>("value"));
    transpCoeff_ = dimensionedScalar("temp", dimensionSet(-1, 0, 2, 0, 0, 1, 0), coeffValue);
}

void transpCoeffsConstValue::calcDiffusionCoeffs()
{
    scalar coeffValue(dict_.get<scalar>("value"));
    transpCoeff_ = dimensionedScalar("temp", dimViscosity, coeffValue);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// He Ion coefficients from LXCat
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
transpCoeffsHeIonLXCat::transpCoeffsHeIonLXCat(const fvMesh& mesh, const dictionary& dict, volScalarField& transpCoeff):
    transportCoeffsHandler(mesh, dict, transpCoeff)
{ }

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void transpCoeffsHeIonLXCat::calcMobilityCoeffs()
{
    // Update EtoN_ field
    EtoN_ = 1.0e21*mag(E_/N_);

    // Coefficients of 4th order polynomial fit from LXCat He Ion mobility data
    scalar a4 = -8.894014448553097e+21;
    scalar a3 = 2.077636242890531e+22;
    scalar a2 = 5.982171417308625e+23;
    scalar a1 = -2.808962742445035e+24;
    scalar a0 = 5.22892285263363e+24;

    transpCoeff_ = dimMu_*(dimN_/N_) * ( a4*Foam::pow(Foam::log(EtoN_/dimEtoN_), 4)
                                       + a3*Foam::pow(Foam::log(EtoN_/dimEtoN_), 3)
                                       + a2*Foam::pow(Foam::log(EtoN_/dimEtoN_), 2)
                                       + a1*Foam::log(EtoN_/dimEtoN_) + a0 );
}

void transpCoeffsHeIonLXCat::calcDiffusionCoeffs()
{
    // Update EtoN_ field
    EtoN_ = 1.0e21*mag(E_/N_);

    // Coefficients of 4th order polynomial fit from LXCat He Ion diffusion data
    scalar a8 = -1.2482977300991656e+16;
    scalar a7 = 1.4550249697805517e+18;
    scalar a6 = -6.9962653253416755e+19;
    scalar a5 = 1.7861620341545414e+21;
    scalar a4 = -2.581844868723774e+22;
    scalar a3 = 2.0689968140002257e+23;
    scalar a2 = -8.477178896864506e+23;
    scalar a1 = 2.1388363733818154e+24;
    scalar a0 = 6.904499705732958e+23;

    transpCoeff_ = dimD_*(dimN_/N_) * ( a8*Foam::pow(EtoN_/dimEtoN_,8) + a7*Foam::pow(EtoN_/dimEtoN_,7)
                                      + a6*Foam::pow(EtoN_/dimEtoN_,6) + a5*Foam::pow(EtoN_/dimEtoN_,5)
                                      + a4*Foam::pow(EtoN_/dimEtoN_,4) + a3*Foam::pow(EtoN_/dimEtoN_,3)
                                      + a2*Foam::pow(EtoN_/dimEtoN_,2) + a1*EtoN_/dimEtoN_ + a0);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Electron coefficients from Pascoa2016
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
transpCoeffsElectronsRelation::transpCoeffsElectronsRelation(const fvMesh& mesh, const dictionary& dict, volScalarField& transpCoeff):
    transportCoeffsHandler(mesh, dict, transpCoeff)
{ }

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void transpCoeffsElectronsRelation::calcMobilityCoeffs()
{
    // Update EtoN_ field
    EtoP_ = mag(E_/P_);

    transpCoeff_ = (dimP_/P_)*dimMu_*(0.21*(24.32*Foam::exp(-(EtoP_/dimEtoP_)/1.057e3) + 19.38*Foam::exp(-(EtoP_/dimEtoP_)/2.3430e4) + 14.45*ones_) +
                                      0.79*(173.1*Foam::exp(-(EtoP_/dimEtoP_)/195.1)   + 36.19*Foam::exp(-(EtoP_/dimEtoP_)/1.2763e4) + 31.73*ones_) );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Diffusion coefficients using the Einstein relation
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
diffusionCoeffsEinsteinRelation::diffusionCoeffsEinsteinRelation(const fvMesh& mesh, const dictionary& dict, volScalarField& transpCoeff, const volScalarField& mobilityCoeff, const scalar chargeNumber):
    transportCoeffsHandler(mesh, dict, transpCoeff),
    mobilityCoeff_(mobilityCoeff),
    chargeNumber_(chargeNumber)
{ }

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void diffusionCoeffsEinsteinRelation::calcDiffusionCoeffs()
{
    scalar tempValue(dict_.get<scalar>("temperature"));

    const dimensionedScalar temp
    (
        "temp",
        dimensionSet(0, 0, 0, 1, 0, 0, 0),
        tempValue
    );

    transpCoeff_ = mobilityCoeff_*constant::physicoChemical::k*temp / (abs(chargeNumber_)*constant::electromagnetic::e) ;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

