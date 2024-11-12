/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 AUTHOR,AFFILIATION
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

#include "parabolicInletFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Foam::scalar Foam::parabolicInletFvPatchVectorField::t() const
// {
//     return db().time().timeOutputValue();
// }


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::parabolicInletFvPatchVectorField::
parabolicInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    center_(0,0,0),
    Um_(1,0,0)


    // fixedValueFvPatchVectorField(p, iF),
    // scalarData_(0),
    // data_(Zero),
    // fieldData_(p.size(), Zero),
    // timeVsData_(),
    // wordData_("wordDefault"),
    // labelData_(-1),
    // boolData_(false)
{
}


Foam::parabolicInletFvPatchVectorField::
parabolicInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    center_(dict.get<vector>("center")),
    Um_(dict.get<vector>("meanVelocity"))

    // fixedValueFvPatchVectorField(p, iF),
    // scalarData_(dict.get<scalar>("scalarData")),
    // data_(dict.get<vector>("data")),
    // fieldData_("fieldData", dict, p.size()),
    // timeVsData_(Function1<vector>::New("timeVsData", dict, &db())),
    // wordData_(dict.getOrDefault<word>("wordName", "wordDefault")),
    // labelData_(-1),
    // boolData_(false)
{


    fixedValueFvPatchVectorField::evaluate();

    /*
    //Initialise with the value entry if evaluation is not possible
    readValueEntry(dict, IOobjectOption::MUST_READ);
    */
}


Foam::parabolicInletFvPatchVectorField::
parabolicInletFvPatchVectorField
(
    const parabolicInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    center_(ptf.center_),
    Um_(ptf.Um_)

    // fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    // scalarData_(ptf.scalarData_),
    // data_(ptf.data_),
    // fieldData_(ptf.fieldData_, mapper),
    // timeVsData_(ptf.timeVsData_.clone()),
    // wordData_(ptf.wordData_),
    // labelData_(-1),
    // boolData_(ptf.boolData_)
{}


Foam::parabolicInletFvPatchVectorField::
parabolicInletFvPatchVectorField
(
    const parabolicInletFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    center_(ptf.center_),
    Um_(ptf.Um_)

    // fixedValueFvPatchVectorField(ptf),
    // scalarData_(ptf.scalarData_),
    // data_(ptf.data_),
    // fieldData_(ptf.fieldData_),
    // timeVsData_(ptf.timeVsData_.clone()),
    // wordData_(ptf.wordData_),
    // labelData_(-1),
    // boolData_(ptf.boolData_)
{}


Foam::parabolicInletFvPatchVectorField::
parabolicInletFvPatchVectorField
(
    const parabolicInletFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    center_(ptf.center_),
    Um_(ptf.Um_)

    // fixedValueFvPatchVectorField(ptf, iF),
    // scalarData_(ptf.scalarData_),
    // data_(ptf.data_),
    // fieldData_(ptf.fieldData_),
    // timeVsData_(ptf.timeVsData_.clone()),
    // wordData_(ptf.wordData_),
    // labelData_(-1),
    // boolData_(ptf.boolData_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::parabolicInletFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
    // fieldData_.autoMap(m);
}


void Foam::parabolicInletFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    // const parabolicInletFvPatchVectorField& tiptf =
    //     refCast<const parabolicInletFvPatchVectorField>(ptf);

    // fieldData_.rmap(tiptf.fieldData_, addr);
}


void Foam::parabolicInletFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // fixedValueFvPatchVectorField::operator==
    // (
    //     data_
    //   + fieldData_
    //   + scalarData_*timeVsData_->value(t())
    // );


    // fixedValueFvPatchVectorField::operator==
    // (
    //     data_
    //   + fieldData_
    //   + scalarData_*timeVsData_->value(t())
    // );


    const vectorField r(patch().Cf() - center_);

    fixedValueFvPatchVectorField::operator==( 2.0*Um_*(1.0-sqr(mag(r)/max(mag(r)))) );
    // fixedValueFvPatchVectorField::operator==( 2.0*Um_*(1.0-sqr(mag(r)/R_)) );

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::parabolicInletFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    os.writeEntry("meanVelocity", Um_);
    os.writeEntry("center", center_);

    // fvPatchVectorField::write(os);
    // os.writeEntry("scalarData", scalarData_);
    // os.writeEntry("data", data_);
    // fieldData_.writeEntry("fieldData", os);
    // timeVsData_->writeData(os);
    // os.writeEntry("wordData", wordData_);
    // fvPatchVectorField::writeValueEntry(os);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        parabolicInletFvPatchVectorField
    );
}

// ************************************************************************* //
