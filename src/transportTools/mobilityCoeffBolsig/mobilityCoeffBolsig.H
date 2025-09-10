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
    Foam::mobilityCoeffBolsig

Description
    This class reads the mobility and diffusion coefficients from a table
    and interpolates them based on the electric field and He concentration.

SourceFiles
    mobilityCoeffBolsig.C

\*---------------------------------------------------------------------------*/

#ifndef mobilityCoeffBolsig_H
#define mobilityCoeffBolsig_H

#include "transportCoeffsBase.H"
#include "interpolation2DTable.H"
#include "tableBounds.H"

namespace Foam
{

class mobilityCoeffBolsig : public transportCoeffsBase
{
protected:
    // Interpolation tables for mobility and diffusion coefficients
    interpolation2DTable<scalar> mobilityTable_;

    // Tabulation variables
    const volScalarField& EN_;
    const volScalarField& He_;

public:
    //- Runtime type information
    TypeName("bolsig");

    // Constructor
    mobilityCoeffBolsig
    (
        const fvMesh& mesh,
        const dictionary& dict,
        volScalarField& transportCoeff
    );

    // Destructor
    // virtual ~transportCoeffsBolsig() {}
    ~mobilityCoeffBolsig() {}

    // Member functions
    //- Calculate mobility coefficients
    void calcMobilityCoeffs();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
