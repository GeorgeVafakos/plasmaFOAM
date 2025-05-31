/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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
    Foam::coupledPotentialSurfaceChargeFvPatchScalarField

Description
    Mixed boundary condition for temperature and radiation heat transfer,
    suitable for multiregion cases. The optional specification of
    thin thermal layer resistances can be made using the entries
    \c thicknessLayers and \c kappaLayers.

    The thermal conductivity \c kappa can be obtained from a variety
    of sources, as explained in the \c temperatureCoupledBase class.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type                compressible::turbulentTemperatureRadCoupledMixed;

        // Optional entries
        Tnbr                <word>;
        qrNbr               <word>;
        qr                  <word>;
        logInterval         <scalar>;
        thermalInertia      <bool>;
        verbose             <bool>;
        prefix              <word>;

        // Conditional entries

            // Option-1
            thicknessLayers <scalarList>;
            kappaLayers     <scalarList>;

            // Option-2
            thicknessLayer  <PatchFunction1<scalar>>;
            kappaLayer      <PatchFunction1<scalar>>;

        // Inherited entries
        kappaMethod         <word>;
        kappa               <word>;
        ...
    }
    \endverbatim

    For example:
    \verbatim
    <patchName>
    {
        type            compressible::turbulentTemperatureRadCoupledMixed;
        Tnbr            T;
        qrNbr           qr; // or none. Name of qr field on neighbour region
        qr              qr; // or none. Name of qr field on this region
        thicknessLayers (0.1 0.2 0.3 0.4);
        kappaLayers     (1 2 3 4);
        thermalInertia  false/true;
        kappaMethod     lookup;
        kappa           kappa;
        value           uniform 300;
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                            | Type  | Reqd | Deflt
      type   | Type name: compressible::turbulentTemperatureRadCoupledMixed <!--
               -->                                      | word  | yes  | -
      Tnbr     | Name of the nbr temperature field      | word  | no   | T
      qrNbr    | Name of the radiative flux in the nbr region | word | no | none
      qr       | Name of the radiative flux in this region | word | no | none
      logInterval | Log-file output frequency [s]          | scalar | no  | -1
      thermalInertia | Flag to add thermal inertia to wall node | bool | no <!--
               -->                                                     | false
      verbose | Flag to enable verbose output with extra fields | bool | no <!--
              -->                                                      | false
      prefix  | Name of output-field prefix (in verbose mode)         <!--
              -->                                       | word | no | multiWorld
      thicknessLayers | List of thicknesses per layer [m] | scalarList <!--
              -->                                             | choice | -
      kappaLayers  | List of thermal conductivites per layer [W/m/K]   <!--
              -->                                | scalarList | choice | -
      thicknessLayer | Single thickness of layer [m]                   <!--
              -->                   | PatchFunction1\<scalar\>| choice | -
      kappaLayer     | Corresponding thermal conductivity [W/m/K]      <!--
              -->                   | PatchFunction1\<scalar\>| choice | -
    \endtable

    The inherited entries are elaborated in:
      - \link mixedFvPatchFields.H \endlink
      - \link temperatureCoupledBase.H \endlink
      - \link mappedPatchFieldBase.H \endlink
      - \link writeFile.H \endlink

Note
  - Needs to be on underlying \c mapped(Wall)FvPatch.

SourceFiles
    coupledPotentialSurfaceChargeFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef coupledPotentialSurfaceChargeFvPatchScalarField_H
#define coupledPotentialSurfaceChargeFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "potentialCoupledBase.H"
#include "mappedPatchFieldBase.H"
#include "scalarField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class coupledPotentialSurfaceChargeFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class coupledPotentialSurfaceChargeFvPatchScalarField 
: 
    public mixedFvPatchScalarField, 
    public potentialCoupledBase,
    public mappedPatchFieldBase<scalar>
{
    // Private data
    
        //- Name of field on the neighbour region
        const word TnbrName_;

        //- Name of the surface charge field
        const word surfaceChargeName_;

    // Private member functions

        //- Calculate coefficients for assembly matrix
        tmp<Field<scalar>> coeffs
        (
            fvMatrix<scalar>& matrix,
            const Field<scalar>&,
            const label
        ) const;


public:

    //- Runtime type information
    TypeName("coupledPotentialSurfaceCharge");


    // Constructors

        //- Construct from patch and internal field
        coupledPotentialSurfaceChargeFvPatchScalarField
        (
            const fvPatch&, 
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        coupledPotentialSurfaceChargeFvPatchScalarField
        (
            const fvPatch&, 
            const DimensionedField<scalar, volMesh>&, 
            const dictionary&
        );

        //- Construct by mapping given
        //- coupledPotentialSurfaceChargeFvPatchScalarField onto a
        //- new patch
        coupledPotentialSurfaceChargeFvPatchScalarField
        (
            const coupledPotentialSurfaceChargeFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        coupledPotentialSurfaceChargeFvPatchScalarField
        (
            const coupledPotentialSurfaceChargeFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        coupledPotentialSurfaceChargeFvPatchScalarField
        (
            const coupledPotentialSurfaceChargeFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return fvPatchField<scalar>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return fvPatchField<scalar>::Clone(*this, iF);
        }


    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<scalar>&,
                const labelList&
            );


    // Evaluation

        //- Given patch potential calculate corresponding epsilon field. Override
        //- potentialCoupledBase::epsilon to includes effect of any
        //- explicit epsilonThickness
        virtual tmp<scalarField> epsilon() const;

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Manipulate matrix
        virtual void manipulateMatrix
        (
            fvMatrix<scalar>& m,
            const label iMatrix,
            const direction cmpt
        );


    // I-O

        //- Write
        virtual void write(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
