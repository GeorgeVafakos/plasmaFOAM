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

#include "ArrheniusRateCoeff.H"
#include "plasmaChemistryModel.H"
#include "volFields.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

void ArrheniusRateCoeff::read(const dictionary& d)
{
    A_ = d.lookupOrDefault<scalar>("A", 1.0);
    B_ = d.lookupOrDefault<scalar>("B", 0.0);
    C_ = d.lookupOrDefault<scalar>("C", 0.0);
    fieldName_ = d.lookupOrDefault<word>("field", "Te");
}

void ArrheniusRateCoeff::calculate(plasmaChemistryModel& chemistry, const label reactionIndex) const
{
    volScalarField& kj = chemistry.k()[reactionIndex];
    const volScalarField& T = chemistry.mesh().lookupObject<volScalarField>("Te");

    kj = dim(kj) * A_*pow(T/dim(T),B_)*exp(-C_/(T/dim(T)));
}




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

