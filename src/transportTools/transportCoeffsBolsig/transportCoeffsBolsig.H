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
    Foam::transportCoeffsBolsig

Description
    This class reads the mobility and diffusion coefficients from a table
    and interpolates them based on the electric field and He concentration.

SourceFiles
    transportCoeffsBolsig.C

\*---------------------------------------------------------------------------*/

#ifndef TRANSPORTCOEFFSBOLSIG_H
#define TRANSPORTCOEFFSBOLSIG_H

#include "transportCoeffsHandler.H"
#include "interpolation2DTable.H"
#include "tableBounds.H"

namespace Foam
{

class transportCoeffsBolsig : public transportCoeffsHandler
{
        interpolation2DTable<scalar> mobilityTable_;
        interpolation2DTable<scalar> diffusivityTable_;

    public:
        //- Runtime type information
        TypeName("transportCoeffsBolsig");

        // Constructor
        transportCoeffsBolsig
        (
            const fvMesh& mesh,
            const dictionary& dict,
            volScalarField& transpCoeff
        );

        // transportCoeffsBolsig
        // (
        //     const fvMesh& mesh,
        //     const dictionary& dict,
        //     volScalarField& electronMobilityField,     // GeometricField<double,...>&
        //     volScalarField& electronDiffusivityField   // if you only pass one, drop this argument
        // );

        // // If solver only passes one field, also add a 3-arg overload:
        // transportCoeffsBolsig
        // (
        //     const fvMesh& mesh,
        //     const dictionary& dict,
        //     volScalarField& electronMobilityField
        // );

        // Destructor
        virtual ~transportCoeffsBolsig() {}
        // ~transportCoeffsBolsig() {}

        // Member functions
        //- Calculate mobility coefficients
        void calcMobilityCoeffs() override;

        //- Calculate diffusion coefficients
        void calcDiffusionCoeffs() override;
};

} // namespace Foam

#endif 