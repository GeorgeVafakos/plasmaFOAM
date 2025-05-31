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

#ifndef transportCoeffsHandler_H
#define transportCoeffsHandler_H

#include "volFields.H"
#include "constants.H"
#include <stdexcept>
#include "IOstreams.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class transportCoeffsHandler Declaration
\*---------------------------------------------------------------------------*/

class transportCoeffsHandler
{
    public:
        //- Runtime type information
        TypeName("transportCoeffsHandler");

        //- Declare run-time constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            transportCoeffsHandler,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict,
                volScalarField& transpCoeff
            ),
            (mesh, dict, transpCoeff)
        );

    protected:
        // fvMesh reference
        const fvMesh& mesh_;

        // Control dictionary reference
        const dictionary& dict_;

        // Transport coefficient reference
        volScalarField& transpCoeff_;

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

        // Ambient pressure
        volScalarField ones_;

    public:
        // Constructor
        transportCoeffsHandler(const fvMesh& mesh, const dictionary& dict, volScalarField& transpCoeff);

        // Destructor
        virtual ~transportCoeffsHandler() {}

        // Member operators
        //- Equality operator
        transportCoeffsHandler operator = (const transportCoeffsHandler& transpCoeffsHandler);
        // voltageHandler& operator=(const voltageHandler& other);

        // Member functions
        //- calculate mobility coefficients
        virtual void calcMobilityCoeffs() { NotImplemented; }

        //- calculate diffusion coefficients
        virtual void calcDiffusionCoeffs() { NotImplemented; }
};


class transpCoeffsConstValue : public transportCoeffsHandler
{
    public:
        // Constructor
        transpCoeffsConstValue(const fvMesh& mesh, const dictionary& dict, volScalarField& transpCoeff);

        // Destructor
        ~transpCoeffsConstValue() {}

        // member functions
        //- calculate mobility coefficients
        void calcMobilityCoeffs();

        //- calculate diffusion coefficients
        void calcDiffusionCoeffs();
};

class transpCoeffsHeIonLXCat : public transportCoeffsHandler
{
    public:
        // Constructor
        transpCoeffsHeIonLXCat(const fvMesh& mesh, const dictionary& dict, volScalarField& transpCoeff);

        // destructor
        ~transpCoeffsHeIonLXCat() {}

        // member functions
        //- calculate mobility coefficients
        void calcMobilityCoeffs();

        //- calculate diffusion coefficients
        void calcDiffusionCoeffs();
};

class transpCoeffsElectronsRelation : public transportCoeffsHandler
{
    public:
        // Constructor
        transpCoeffsElectronsRelation(const fvMesh& mesh, const dictionary& dict, volScalarField& transpCoeff);

        // Destructor
        ~transpCoeffsElectronsRelation() {}

        // member functions
        //- calculate mobility coefficients
        void calcMobilityCoeffs();

        //- calculate diffusion coefficients
        // void calcDiffusionCoeffs();
};

class diffusionCoeffsEinsteinRelation : public transportCoeffsHandler
{
    protected:
        // Mobility of species
        const volScalarField& mobilityCoeff_;

        // Charge number of species
        const scalar chargeNumber_;

    public:
        // Constructor
        diffusionCoeffsEinsteinRelation(const fvMesh& mesh, const dictionary& dict, volScalarField& transpCoeff, const volScalarField& mobilityCoeff, const scalar chargeNumber);

        // Destructor
        ~diffusionCoeffsEinsteinRelation() {}

        // member functions
        //- calculate mobility coefficients
        // void calcMobilityCoeffs();

        //- calculate diffusion coefficients
        void calcDiffusionCoeffs();
};





// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
