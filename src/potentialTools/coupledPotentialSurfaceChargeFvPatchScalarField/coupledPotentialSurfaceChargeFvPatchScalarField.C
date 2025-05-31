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

#include "coupledPotentialSurfaceChargeFvPatchScalarField.H"
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
coupledPotentialSurfaceChargeFvPatchScalarField::
coupledPotentialSurfaceChargeFvPatchScalarField
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
//- coupledPotentialSurfaceChargeFvPatchScalarField onto a
//- new patch
coupledPotentialSurfaceChargeFvPatchScalarField::
coupledPotentialSurfaceChargeFvPatchScalarField
(
    const coupledPotentialSurfaceChargeFvPatchScalarField& ptf,
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
coupledPotentialSurfaceChargeFvPatchScalarField::
coupledPotentialSurfaceChargeFvPatchScalarField
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
    TnbrName_(dict.get<word>("Tnbr")),
    surfaceChargeName_(dict.getOrDefault<word>("surfaceCharge", "surfaceCharge"))
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
coupledPotentialSurfaceChargeFvPatchScalarField::
coupledPotentialSurfaceChargeFvPatchScalarField
(
    const coupledPotentialSurfaceChargeFvPatchScalarField& wtcsf,
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
coupledPotentialSurfaceChargeFvPatchScalarField::
coupledPotentialSurfaceChargeFvPatchScalarField
(
    const coupledPotentialSurfaceChargeFvPatchScalarField& wtcsf
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

void coupledPotentialSurfaceChargeFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& mapper
)
{
    mixedFvPatchScalarField::autoMap(mapper);
    potentialCoupledBase::autoMap(mapper);
    //mappedPatchFieldBase<scalar>::autoMap(mapper);
}

void coupledPotentialSurfaceChargeFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const coupledPotentialSurfaceChargeFvPatchScalarField& tiptf =
        refCast<const coupledPotentialSurfaceChargeFvPatchScalarField>(ptf);

    potentialCoupledBase::rmap(tiptf, addr);
    //mappedPatchFieldBase<scalar>::rmap(ptf, addr);
}

tmp<scalarField> coupledPotentialSurfaceChargeFvPatchScalarField::epsilon() const
{
    return potentialCoupledBase::epsilon();
}

void coupledPotentialSurfaceChargeFvPatchScalarField::updateCoeffs()
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
    const tmp<scalarField> myEpsilonDelta = epsilonVp * patch().deltaCoeffs(); // Print deltaCoeffs to see what they are

    scalarField nbrIntFld;
    scalarField nbrEpsilonDelta;
    if (mpp.sameWorld())
    {
        // Same world so lookup
        const auto& nbrMesh = refCast<const fvMesh>(mpp.sampleMesh());
        const label nbrPatchID = mpp.samplePolyPatch().index();
        const auto& nbrPatch = nbrMesh.boundary()[nbrPatchID];

        const coupledPotentialSurfaceChargeFvPatchScalarField& nbrField =
            refCast<const coupledPotentialSurfaceChargeFvPatchScalarField>
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

    // Add surface charge contribution
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const volScalarField& surfaceCharge =
        mesh.lookupObject<volScalarField>(surfaceChargeName_);

    const scalarField& surfaceChargePatch = surfaceCharge.boundaryField()[patch().index()];

    this->refValue() = nbrIntFld;
    this->refGrad() = surfaceChargePatch/epsilonVp;
    this->valueFraction() = nbrEpsilonDelta / (nbrEpsilonDelta + myEpsilonDelta());

    mixedFvPatchScalarField::updateCoeffs();
}

void coupledPotentialSurfaceChargeFvPatchScalarField::manipulateMatrix
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

void coupledPotentialSurfaceChargeFvPatchScalarField::write(Ostream& os) const
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
    coupledPotentialSurfaceChargeFvPatchScalarField
);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //












