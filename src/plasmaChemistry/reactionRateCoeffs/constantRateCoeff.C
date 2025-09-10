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

#include "constantRateCoeff.H"
#include "plasmaChemistryModel.H"
#include "volFields.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(constantRateCoeff, 0);
addToRunTimeSelectionTable(reactionRateCoeffsBase, constantRateCoeff, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
constantRateCoeff::constantRateCoeff(const dictionary& dict, const word& name, plasmaChemistryModel& chemistry)
:
    constantValue_(dict.lookupOrDefault<scalar>("value", 1e-10))
{ 
    isConstant_ = true;  // explicitly set for none value
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
// void constantRateCoeff::read(const dictionary& d)
// {
//     // constantValue_ = d.lookupOrDefault<scalar>("value", 1e-10);
// }

void constantRateCoeff::calculate(plasmaChemistryModel& chemistry, const label reactionIndex) const
{
    volScalarField& k = chemistry.k()[reactionIndex];
    k = dim(k) * constantValue_;
    // k.correctBoundaryConditions(); // Optional
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

