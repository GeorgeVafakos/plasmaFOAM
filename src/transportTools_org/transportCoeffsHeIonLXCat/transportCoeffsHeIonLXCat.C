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

#include "transportCoeffsHeIonLXCat.H" 
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(transportCoeffsHeIonLXCat, 0);
addToRunTimeSelectionTable(transportCoeffsBase, transportCoeffsHeIonLXCat, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
transportCoeffsHeIonLXCat::transportCoeffsHeIonLXCat(const fvMesh& mesh, const dictionary& dict, volScalarField& transportCoeff):
    transportCoeffsBase(mesh, dict, transportCoeff),
    fieldName_(dict_.lookupOrDefault<word>("field", "EN")),
    EN_(mesh_.lookupObjectRef<volScalarField>(fieldName_))
{ }

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void transportCoeffsHeIonLXCat::calcMobilityCoeffs()
{
    // Coefficients of 4th order polynomial fit from LXCat He Ion mobility data
    scalar a4 = -8.894014448553097e+21;
    scalar a3 = 2.077636242890531e+22;
    scalar a2 = 5.982171417308625e+23;
    scalar a1 = -2.808962742445035e+24;
    scalar a0 = 5.22892285263363e+24;

    transportCoeff_ = dimMu_*( a4*Foam::pow(Foam::log(nonDim(EN_)), 4)
                             + a3*Foam::pow(Foam::log(nonDim(EN_)), 3)
                             + a2*Foam::pow(Foam::log(nonDim(EN_)), 2)
                             + a1*Foam::log(nonDim(EN_)) + a0 );

    transportCoeff_.clamp_min(2.05e+24);
    transportCoeff_.clamp_max(2.789e+25);

    transportCoeff_ /= nonDim(N_);
}

void transportCoeffsHeIonLXCat::calcDiffusionCoeffs()
{
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

    transportCoeff_ = dimD_*( a8*Foam::pow(nonDim(EN_),8) + a7*Foam::pow(nonDim(EN_),7)
                            + a6*Foam::pow(nonDim(EN_),6) + a5*Foam::pow(nonDim(EN_),5)
                            + a4*Foam::pow(nonDim(EN_),4) + a3*Foam::pow(nonDim(EN_),3)
                            + a2*Foam::pow(nonDim(EN_),2) + a1*nonDim(EN_) + a0);

    transportCoeff_.clamp_min(6.42e+23);
    transportCoeff_.clamp_max(1.1e+25);

    transportCoeff_ /= nonDim(N_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

