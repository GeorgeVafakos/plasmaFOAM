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

#include "mobilityCoeffUniformValue.H" 
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(mobilityCoeffUniformValue, 0);
addToRunTimeSelectionTable(transportCoeffsBase, mobilityCoeffUniformValue, mobility);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
mobilityCoeffUniformValue::mobilityCoeffUniformValue(const fvMesh& mesh, const dictionary& dict, volScalarField& transportCoeff):
    transportCoeffsBase(mesh, dict, transportCoeff)
{
    isConstant_ = true;  // explicitly set for uniform value
    calcMobilityCoeffs();    // calculate once upon construction
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void mobilityCoeffUniformValue::calcMobilityCoeffs()
{
    scalar coeffValue(dict_.get<scalar>("value"));
    transportCoeff_ = dimensionedScalar("temp", dimensionSet(-1, 0, 2, 0, 0, 1, 0), coeffValue);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

