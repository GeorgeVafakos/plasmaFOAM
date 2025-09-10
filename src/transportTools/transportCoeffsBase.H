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

Class
    Foam::transportCoeffsHandler

Description
    This is a class that provides the calculation of the transport
    coefficients (mobility and diffusion) for the ions.

SourceFiles
    transportCoeffsHandler.C

\*---------------------------------------------------------------------------*/

#ifndef transportCoeffsBase_H
#define transportCoeffsBase_H

#include "volFields.H"
#include "constants.H"
#include "IOstreams.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class transportCoeffsBase Declaration
\*---------------------------------------------------------------------------*/

class transportCoeffsBase
{
protected:
    // fvMesh reference
    const fvMesh& mesh_;

    // Control dictionary reference
    const dictionary& dict_;

    // Transport coefficient reference
    volScalarField& transportCoeff_;

    // Electric field
    const volVectorField& E_;

    // Number density
    const volScalarField& N_;

    // Ambient pressure
    const volScalarField& P_;

    // Reduced electric field (Td)
    volScalarField EtoN_;

    // Eelectric field over pressure
    volScalarField EtoP_;

    // Unit dimentioned scalars
    const dimensionedScalar dimEtoN_;
    const dimensionedScalar dimEtoP_;
    const dimensionedScalar dimMu_;
    const dimensionedScalar dimD_;
    const dimensionedScalar dimN_;
    const dimensionedScalar dimP_;

    // Flag indicating constant coefficient throughout the simulation
    bool isConstant_;

    //- inline - return a non-dimensioned volScalarField field
    volScalarField nonDim(const volScalarField& field) const
    {
        return field/dimensionedScalar("field", field.dimensions(), scalar(1.0));
    }

public:
    //- Runtime type information
    TypeName("transportCoeffsBase");

    // Declare run-time selection table for mobility coefficients
    declareRunTimeSelectionTable
    (
        autoPtr,
        transportCoeffsBase,
        mobility,
        (
            const fvMesh& mesh,
            const dictionary& dict,
            volScalarField& transportCoeff
        ),
        (mesh, dict, transportCoeff)
    );

    // Declare run-time selection table for diffusion coefficients
    declareRunTimeSelectionTable
    (
        autoPtr,
        transportCoeffsBase,
        diffusion,
        (
            const fvMesh& mesh,
            const dictionary& dict,
            volScalarField& transportCoeff
        ),
        (mesh, dict, transportCoeff)
    );

    // Selectors
    //- Return a reference to the selected coefficient 
    static autoPtr<transportCoeffsBase> NewMobility
    (
        const fvMesh& mesh,
        const dictionary& dict,
        volScalarField& transportCoeff
    );

    static autoPtr<transportCoeffsBase> NewDiffusion
    (
        const fvMesh& mesh,
        const dictionary& dict,
        volScalarField& transportCoeff
    );

    // Constructor
    transportCoeffsBase(const fvMesh& mesh, const dictionary& dict, volScalarField& transportCoeff);

    // Destructor
    virtual ~transportCoeffsBase() {}

    // Member operators
    //- Equality operator
    transportCoeffsBase operator = (const transportCoeffsBase& transpCoeffsBase);
    // voltageHandler& operator=(const voltageHandler& other);

    // Member functions
    //- calculate mobility coefficients
    virtual void calcMobilityCoeffs() { NotImplemented; }

    //- calculate diffusion coefficients
    virtual void calcDiffusionCoeffs() { NotImplemented; }

    //- inline - return true to ignore calculation of coefficients
    bool isConstant() const 
    {
        return isConstant_;
    }
};






// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
