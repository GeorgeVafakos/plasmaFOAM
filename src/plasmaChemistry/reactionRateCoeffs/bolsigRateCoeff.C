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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
bolsigRateCoeff::bolsigRateCoeff(const dictionary& dict, const word& name, plasmaChemistryModel& chemistry)
:
    bolsigProperties_
    (
        IOobject
        (
            "bolsigProperties",
            chemistry.mesh().time().constant(),
            chemistry.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    )
{
    scalarList ENList_
    (
        bolsigProperties_.lookup("EN")
    );

    scalarList HeList_
    (
        bolsigProperties_.lookup("He")
    );

    List<List<scalar>> rateCoefficientsList_
    (
        bolsigProperties_.lookup("rateCoefficients")
    );
    
    // Initialise pointer list size
    ratesTable_.setSize(rateCoefficientsList_[0].size());

    forAll(rateCoefficientsList_[0], r)
    {
        // Create lookup table for each reaction
        List<Tuple2<scalar, List<Tuple2<scalar, scalar>>>> rTable(HeList_.size());

        forAll(HeList_, i)
        {
            List<Tuple2<scalar, scalar>> rRow(ENList_.size());

            forAll(ENList_, j)
            {
                label index = i * ENList_.size() + j;
                rRow[j]  = Tuple2<scalar, scalar>(ENList_[j], rateCoefficientsList_[index][r]);
            }
            rTable[i]  = Tuple2<scalar, List<Tuple2<scalar, scalar>>>(HeList_[i], rRow);
        }

        ratesTable_.set
        (
            r,
            new interpolation2DTable<scalar>
            (
                rTable,
                bounds::normalBounding(bounds::normalBounding::CLAMP),
                fileName("ratesTable")
            )
        );
    }
}


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
    const volScalarField& EN_ = chemistry.mesh().lookupObjectRef<volScalarField>("E_N");
    const volScalarField& HeF_ = chemistry.mesh().lookupObjectRef<volScalarField>("HeF");


    volScalarField& kj = chemistry.k()[j];

    // kj = dim(kj) * A_*pow(T/dim(T),B_)*exp(C_/(T/dim(T)));


    

    forAll(kj, celli)
    {
        kj[celli] = ratesTable_[j](HeF_[celli], EN_[celli]);
    }



}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

