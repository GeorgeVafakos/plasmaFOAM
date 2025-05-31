/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "potentialCoupledBase.H"
#include "volFields.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Default constructor
potentialCoupledBase::potentialCoupledBase
(
    const fvPatch& patch, 
    const word& epsilonName
) 
: 
    patch_(patch), 
    epsilonName_(epsilonName)
{}  

//- Construct from patch and dictionary
potentialCoupledBase::potentialCoupledBase
(
    const fvPatch& patch, 
    const dictionary& dict
) 
: 
    patch_(patch), 
    epsilonName_(dict.getOrDefault<word>("epsilon", "epsilon"))
{}

//- Construct from patch and potentialCoupledBase
potentialCoupledBase::potentialCoupledBase
(
    const fvPatch& patch, 
    const potentialCoupledBase& base
) 
: 
    patch_(patch), 
    epsilonName_(base.epsilonName_)
{}

//- Copy construct
potentialCoupledBase::potentialCoupledBase
(
    const potentialCoupledBase& base
) 
:
    potentialCoupledBase(base.patch_, base)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void potentialCoupledBase::autoMap
(
    const fvPatchFieldMapper& mapper
)
{}

void potentialCoupledBase::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{}

tmp<scalarField> potentialCoupledBase::epsilon() const
{
    const fvMesh& mesh = patch_.boundaryMesh().mesh();
    // const label patchi = patch_.index();

    {
        const auto* ptr = mesh.cfindObject<volScalarField>(epsilonName_);

        if (ptr)
        {
            return patch_.patchField<volScalarField>(*ptr);
        }
    }

    {
        const auto* ptr =
            mesh.cfindObject<volSymmTensorField>(epsilonName_);

        if (ptr)
        {
            const symmTensorField& wallValues =
                patch_.patchField<volSymmTensorField>(*ptr);

            const vectorField n(patch_.nf());

            return n & wallValues & n;
        }
    }

    // Error handling if field not found
    FatalErrorInFunction 
        << "Did not find field '" << epsilonName_
        << "' on mesh " << mesh.name() 
        << " patch " << patch_.name()
        << "Please set 'kappa' to the name of"
           " a volScalar or volSymmTensor field"
           ", or use another method" << nl
        << exit(FatalError);

    return scalarField();
}

void potentialCoupledBase::write(Ostream& os) const
{
    os.writeEntry("epsilon", epsilonName_);
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

