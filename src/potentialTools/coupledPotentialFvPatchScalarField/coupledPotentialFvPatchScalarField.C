/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "coupledPotentialFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"

#include "IOField.H"
#include "mappedPatchFieldBase.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from patch and internal field
coupledPotentialFvPatchScalarField::
coupledPotentialFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    potentialCoupledBase(patch()),  // default method (fluidThermo)
    mappedPatchFieldBase<scalar>
    (
        mappedPatchFieldBase<scalar>::mapper(p, iF),
        *this
    ),
    TnbrName_("undefined-Tnbr")
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 1.0;
}

//- Construct by mapping given
//- coupledPotentialFvPatchScalarField onto a
//- new patch
coupledPotentialFvPatchScalarField::
coupledPotentialFvPatchScalarField
(
    const coupledPotentialFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    potentialCoupledBase(patch(), ptf),
    mappedPatchFieldBase<scalar>
    (
        mappedPatchFieldBase<scalar>::mapper(p, iF),
        *this,
        ptf
    ),
    TnbrName_(ptf.TnbrName_)
{}

//- Construct from patch, internal field and dictionary
coupledPotentialFvPatchScalarField::
coupledPotentialFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    potentialCoupledBase(p, dict),
    mappedPatchFieldBase<scalar>
    (
        mappedPatchFieldBase<scalar>::mapper(p, iF),
        *this,
        dict
    ),
    TnbrName_(dict.get<word>("Tnbr"))
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorInFunction
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }

    this->readValueEntry(dict, IOobjectOption::MUST_READ);

    if (this->readMixedEntries(dict))
    {
        // Full restart
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = Zero;
        valueFraction() = 1.0;
    }
}


//- Construct as copy setting internal field reference
coupledPotentialFvPatchScalarField::
coupledPotentialFvPatchScalarField
(
    const coupledPotentialFvPatchScalarField& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(wtcsf, iF),
    potentialCoupledBase(patch(), wtcsf),
    mappedPatchFieldBase<scalar>
    (
        mappedPatchFieldBase<scalar>::mapper(patch(), iF),
        *this,
        wtcsf
    ),
    TnbrName_(wtcsf.TnbrName_)
{}

//- Construct as copy setting internal field reference
coupledPotentialFvPatchScalarField::
coupledPotentialFvPatchScalarField
(
    const coupledPotentialFvPatchScalarField& wtcsf
)
:
    mixedFvPatchScalarField(wtcsf),
    potentialCoupledBase(patch(), wtcsf),
    mappedPatchFieldBase<scalar>
    (
        mappedPatchFieldBase<scalar>::mapper(patch(), wtcsf.internalField()),
        *this,
        wtcsf
    ),
    TnbrName_(wtcsf.TnbrName_)
{}






// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void coupledPotentialFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& mapper
)
{
    mixedFvPatchScalarField::autoMap(mapper);
    potentialCoupledBase::autoMap(mapper);
    //mappedPatchFieldBase<scalar>::autoMap(mapper);
}

void coupledPotentialFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const coupledPotentialFvPatchScalarField& tiptf =
        refCast<const coupledPotentialFvPatchScalarField>(ptf);

    potentialCoupledBase::rmap(tiptf, addr);
    //mappedPatchFieldBase<scalar>::rmap(ptf, addr);
}

tmp<scalarField> coupledPotentialFvPatchScalarField::epsilon() const
{
    return potentialCoupledBase::epsilon();
}

void coupledPotentialFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        mappedPatchFieldBase<scalar>::mapper
        (
            patch(),
            this->internalField()
        );

    // const scalarField& Vp = *this;
    const scalarField epsilonVp(epsilon());
    const tmp<scalarField> myEpsilonDelta = epsilonVp * patch().deltaCoeffs();

    scalarField nbrIntFld;
    scalarField nbrEpsilonDelta;
    if (mpp.sameWorld())
    {
        // Same world so lookup
        const auto& nbrMesh = refCast<const fvMesh>(mpp.sampleMesh());
        const label nbrPatchID = mpp.samplePolyPatch().index();
        const auto& nbrPatch = nbrMesh.boundary()[nbrPatchID];

        const coupledPotentialFvPatchScalarField& nbrField =
            refCast<const coupledPotentialFvPatchScalarField>
            (
                nbrPatch.lookupPatchField<volScalarField>(TnbrName_)
            );


        // Swap to obtain full local values of neighbour epsilon*deltaCoeffs
        nbrIntFld = nbrField.patchInternalField();
        nbrEpsilonDelta = nbrField.epsilon() * nbrPatch.deltaCoeffs();
    }
    else
    {
        // Different world so use my region,patch. Distribution below will
        // do the reordering.
        nbrIntFld = patchInternalField();
        nbrEpsilonDelta = myEpsilonDelta.ref();
    }
    distribute(this->internalField().name() + "_value", nbrIntFld);
    distribute(this->internalField().name() + "_weights", nbrEpsilonDelta);

    this->refValue() = nbrIntFld;
    this->refGrad() = Zero;
    this->valueFraction() = nbrEpsilonDelta / (nbrEpsilonDelta + myEpsilonDelta());

    mixedFvPatchScalarField::updateCoeffs();

    // UPstream::msgType(oldTag);  // Restore tag
}

void coupledPotentialFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& m,
    const label iMatrix,
    const direction cmpt
)
{
    FatalErrorInFunction
        << "This BC does not support matrix manipulation"
        << abort(FatalError);
}

void coupledPotentialFvPatchScalarField::write(Ostream& os) const
{
    mixedFvPatchField<scalar>::write(os);
    os.writeEntry("Tnbr", TnbrName_);
    potentialCoupledBase::write(os);
    mappedPatchFieldBase<scalar>::write(os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    coupledPotentialFvPatchScalarField
);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //












