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
    Foam::voltageHandler

Description
    This is a class that provides the calculation of the electric 
    potential that corresponds to the external electric field, in
    the current time step.

    There are currently two modes of operation:
        - AC
        - Pulsed
        - DC

SourceFiles
    voltageHandler.C

\*---------------------------------------------------------------------------*/

#ifndef voltageHandler_H
#define voltageHandler_H

#include "fvMesh.H"
#include "volFields.H"
#include <stdexcept>
#include "IOstreams.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class voltageHandler Declaration
\*---------------------------------------------------------------------------*/

class voltageHandler
{
    protected:
        // fvMesh reference
        const fvMesh& mesh_;

        // Dictionary reference
        const dictionary& dict_;

        // Potential of the external electric field
        volScalarField& voltExt_;

        // Potential of the external electric field amplitude
        const volScalarField& voltExtAmp_;

    public:
        // Constructor
        voltageHandler(const fvMesh& mesh, const dictionary& dict, volScalarField& voltExt, const volScalarField& voltExtAmp);

        // Destructor
        virtual ~voltageHandler() {}

        // Member operators
        //- Equality operator
        voltageHandler operator = (const voltageHandler& voltHandler);
        // voltageHandler& operator=(const voltageHandler& other);

        // Member functions
        //- calculate voltage field
        virtual volScalarField calcVoltage() { return voltExt_; }
};


class voltageAC : public voltageHandler
{
    public:
        // Constructor
        voltageAC(const fvMesh& mesh, const dictionary& dict, volScalarField& voltExt, const volScalarField& voltExtAmp);

        // destructor
        ~voltageAC() {}

        // member functions
        //- calculate voltage field
        volScalarField calcVoltage();

};


class voltagePulsed : public voltageHandler
{
    public:
        // constructor
        voltagePulsed(const fvMesh& mesh, const dictionary& dict, volScalarField& voltExt, const volScalarField& voltExtAmp);

        // destructor
        ~voltagePulsed() {}

        // Member Functions
        //- calculate voltage field
        volScalarField calcVoltage();
};


class voltageDC : public voltageHandler
{
    public:

        // constructor
        voltageDC(const fvMesh& mesh, const dictionary& dict, volScalarField& voltExt, const volScalarField& voltExtAmp);

        // destructor
        ~voltageDC() {}

        // Member Functions
        //- calculate voltage field
        volScalarField calcVoltage();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
