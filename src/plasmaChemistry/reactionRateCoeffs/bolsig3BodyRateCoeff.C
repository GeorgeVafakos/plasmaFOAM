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

#include "bolsig3BodyRateCoeff.H"
#include "addToRunTimeSelectionTable.H"
#include "plasmaChemistryModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(bolsig3BodyRateCoeff, 0);
addToRunTimeSelectionTable(reactionRateCoeffsBase, bolsig3BodyRateCoeff, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
bolsig3BodyRateCoeff::bolsig3BodyRateCoeff
(
    const dictionary& dict,
    const word& name,
    plasmaChemistryModel& chemistry
)
:
    Foam::bolsigRateCoeff(dict, name, chemistry),
    NGas_(chemistry.mesh().lookupObjectRef<volScalarField>("NGas"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void bolsig3BodyRateCoeff::calculate(plasmaChemistryModel& chemistry, const label j) const
{
    if (!HePtr_)
    {
        HePtr_ = &chemistry.mesh().lookupObjectRef<volScalarField>("He");
    }

    const volScalarField& He = *HePtr_;
    volScalarField& kj = chemistry.k()[j];

    forAll(kj, celli)
    {
        scalar twoBodyRate = ratesTable_[j](He[celli], EN_[celli]);
        kj[celli] = twoBodyRate * NGas_[celli];
    }
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

