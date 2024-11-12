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

Class
    Foam::parabolicInletFvPatchVectorField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a parabolicInlet condition,
    calculated as:

        \f[
            Q = Q_{0} + Q_{p} + s*Q_{t}
        \f]

    where
    \vartable
        s      | single scalar value [units]
        Q_{0}  | single vector value [units]
        Q_{p}  | vector field across patch [units]
        Q_{t}  | vector function of time [units]
    \endtable

Usage
    \table
        Property   | Description                  | Req'd | Default
        scalarData | single scalar value          | yes |
        data       | single vector value          | yes |
        fieldData  | vector field across patch    | yes |
        timeVsData | vector function of time      | yes |
        wordData   | word, eg name of data object | no  | wordDefault
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type       parabolicInlet;
        scalarData -1;
        data       (1 0 0);
        fieldData  uniform (3 0 0);
        timeVsData table (
                             (0 (0 0 0))
                             (1 (2 0 0))
                          );
        wordName   anotherName;
        value      uniform (4 0 0); // optional initial value
    }
    \endverbatim

SourceFiles
    parabolicInletFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef parabolicInletFvPatchVectorField_H
#define parabolicInletFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class parabolicInletFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class parabolicInletFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Single valued scalar quantity, e.g. a coefficient
        // scalar R_;

        //- Single valued Type quantity, e.g. reference pressure pRefValue_
        //  Other options include vector, tensor
        vector center_;
        vector Um_;

        

    //     //- Single valued Type quantity, e.g. reference pressure pRefValue_
    //     //  Other options include vector, tensor
    //     vector data_;

    //     //- Field of Types, typically defined across patch faces
    //     //  e.g. total pressure p0_.  Other options include vectorField
    //     vectorField fieldData_;

    //     //- Type specified as a function of time for time-varying BCs
    //     autoPtr<Function1<vector>> timeVsData_;

    //     //- Word entry, e.g. pName_ for name of the pressure field on database
    //     word wordData_;

    //     //- Label, e.g. patch index, current time index
    //     label labelData_;

    //     //- Boolean for true/false, e.g. specify if flow rate is volumetric_
    //     bool boolData_;


    // Private Member Functions

    //     //- Return current time
    //     scalar t() const;


public:

    //- Runtime type information
    TypeName("parabolicInlet");


    // Constructors

        //- Construct from patch and internal field
        parabolicInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        parabolicInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        parabolicInletFvPatchVectorField
        (
            const parabolicInletFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy construct
        parabolicInletFvPatchVectorField
        (
            const parabolicInletFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new parabolicInletFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        parabolicInletFvPatchVectorField
        (
            const parabolicInletFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new parabolicInletFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
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
                const fvPatchVectorField&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
