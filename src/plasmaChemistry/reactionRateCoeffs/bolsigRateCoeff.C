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

#include "bolsigRateCoeff.H"
#include "plasmaChemistryModel.H"
#include "volFields.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(bolsigRateCoeff, 0);
    addToRunTimeSelectionTable(reactionRateCoeffsBase, bolsigRateCoeff, dictionary);


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void bolsigRateCoeff::read(const dictionary& dict)
{
    A_ = dict.lookupOrDefault<scalar>("A", 1.0);
    B_ = dict.lookupOrDefault<scalar>("B", 0.0);
    C_ = dict.lookupOrDefault<scalar>("C", 0.0);
    fieldName_ = dict.lookupOrDefault<word>("field", "Te");
}

void bolsigRateCoeff::calculate(plasmaChemistryModel& chemistry, const label j) const
{
    volScalarField& kj = chemistry.k()[j];
    const volScalarField& T = chemistry.mesh().lookupObject<volScalarField>("Te");

    kj = dim(kj) * A_*pow(T/dim(T),B_)*exp(C_/(T/dim(T)));
}




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

