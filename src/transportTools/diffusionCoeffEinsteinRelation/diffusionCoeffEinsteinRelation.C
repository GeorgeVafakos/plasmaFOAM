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

#include "diffusionCoeffEinsteinRelation.H" 
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(diffusionCoeffEinsteinRelation, 0);
addToRunTimeSelectionTable(transportCoeffsBase, diffusionCoeffEinsteinRelation, diffusion);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
diffusionCoeffEinsteinRelation::diffusionCoeffEinsteinRelation(const fvMesh& mesh, const dictionary& dict, volScalarField& transportCoeff):
    transportCoeffsBase(mesh, dict, transportCoeff),
    fieldName_(dict_.lookupOrDefault<word>("T", "TEle")),
    chargeNumber_(dict.parent().get<int>("chargeNumber")),
    T_(mesh_.lookupObjectRef<volScalarField>(fieldName_)),
    mobilityCoeff_(mesh_.lookupObjectRef<volScalarField>("mobilityCoeffSpecies_"+dict_.parent().dictName())),
    EinsteinFactor_(constant::physicoChemical::k / (abs(chargeNumber_)*constant::electromagnetic::e))
{ }

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
// void diffusionCoeffEinsteinRelation::calcMobilityCoeffs()
// {
//     FatalErrorInFunction
//         << "The diffusionEinsteinRelation can not be used "
//         << "for the mobility of species: " << dict_.parent().dictName()
//         << abort(FatalError);
// }

void diffusionCoeffEinsteinRelation::calcDiffusionCoeffs()
{
    transportCoeff_ = mobilityCoeff_ * (EinsteinFactor_ * T_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

