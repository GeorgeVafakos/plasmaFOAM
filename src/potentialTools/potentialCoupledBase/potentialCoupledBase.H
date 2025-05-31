/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019,2021 OpenCFD Ltd.
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

Class
    Foam::potentialCoupledBase

Description
    Common functions used in electric potential coupled boundaries.

    The electric permittivity \c kappa may be obtained by the following methods:
      - 'lookup' : lookup volScalarField (or volSymmTensorField) with name
        defined by 'kappa'
      - 'fluidThermo' : use fluidThermo and default
        compressible::turbulenceModel to calculate kappa
      - 'solidThermo' : use solidThermo kappa()
      - 'directionalSolidThermo': uses look up for volSymmTensorField for
        transformed kappa vector. Field name definable in 'alphaAni',
        named 'Anialpha' in solid solver by default
      - 'function' : kappa, alpha directly specified as Function1
      - 'phaseSystem' : used for multiphase thermos

    \par Keywords provided by this class:
      \table
        Property     | Description                          | Required | Default
        kappaMethod  | Thermal conductivity method          | yes |
        kappa        | Name of thermal conductivity field   | partly |
        alpha        | Name of thermal diffusivity field    | partly |
        alphaAni     | Name of non-isotropic alpha          | partly |
        kappaValue   | Function1 supplying kappa            | partly |
        alphaValue   | Function1 supplying alpha            | partly |
      \endtable

Usage
    \verbatim
    nonIsotropicWall
    {
        ...
        kappaMethod     directionalSolidThermo;
        alphaAni        Anialpha;
        ...
    }
    \endverbatim

    \verbatim
    specifiedWall
    {
        ...
        kappaMethod     function;
        kappaFunction   constant 1.0;
        alphaFunction   constant 100.0;
        ...
    }
    \endverbatim

SourceFiles
    temperatureCoupledBase.C

\*---------------------------------------------------------------------------*/

#ifndef potentialCoupledBase_H
#define potentialCoupledBase_H

#include "scalarField.H"
#include "fvPatch.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class potentialCoupledBase Declaration
\*---------------------------------------------------------------------------*/

class potentialCoupledBase
{
protected:

    // Protected Data

        //- Underlying patch (e.g. boundary patch reference)
        const fvPatch& patch_;

        //- Name of electric permittivity field (if looked up from database)
        const word epsilonName_;

public:

    // Constructors

        //- Default constructor
        explicit potentialCoupledBase
        (
          const fvPatch& patch, 
          const word& epsilonName = "epsilon"
        );

        //- Construct from patch and dictionary
        potentialCoupledBase
        (
          const fvPatch& patch, 
          const dictionary& dict
        );

        //- Construct from patch and potentialCoupledBase
        potentialCoupledBase
        (
          const fvPatch& patch, 
          const potentialCoupledBase& base
        );

        //- Copy construct
        potentialCoupledBase(const potentialCoupledBase& base);


    // Destructor
    virtual ~potentialCoupledBase() = default;

    
    // Member Functions

        //- Name of electric permittivity field
        const word& epsilonName() const noexcept 
        {
          return epsilonName_; 
        }

        //- Map (and resize as needed) from self given a mapping object
        virtual void autoMap
        (
          const fvPatchFieldMapper&
        ) = 0;

        //- Reverse map the given fvPatchField onto this fvPatchField
        virtual void rmap
        (
          const fvPatchField<scalar>&, 
          const labelList&
        ) = 0;

        //- Calculate corresponding electric permittivity field
        virtual tmp<scalarField> epsilon() const;

        //- Write
        void write(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
