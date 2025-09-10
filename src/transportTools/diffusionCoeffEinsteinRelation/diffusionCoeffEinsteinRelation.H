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
    Foam::transportCoeffsUniformValue

Description
    This is a class that provides the calculation of the transport
    coefficients (mobility and diffusion) for the ions.

SourceFiles
    transportCoeffsUniformValue.C

\*---------------------------------------------------------------------------*/

#ifndef transportCoeffsUniformValue_H
#define transportCoeffsUniformValue_H

#include "transportCoeffsBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class transportCoeffsUniformValue Declaration
\*---------------------------------------------------------------------------*/

class diffusionCoeffEinsteinRelation : public transportCoeffsBase
{
protected:
    // Name of the temperature field
    word fieldName_;

    // Charge number
    int chargeNumber_;

    // Tempreture field
    volScalarField& T_;

    // Mobility coefficient of species
    volScalarField& mobilityCoeff_;

    // Einstein factor (precalculate coefficients to make calculations faster)
    dimensionedScalar EinsteinFactor_;

public:
    TypeName("diffusionEinsteinRelation");

    // Constructor
    diffusionCoeffEinsteinRelation
    (
        const fvMesh& mesh, 
        const dictionary& dict, 
        volScalarField& transportCoeff
    );

    // Destructor
    ~diffusionCoeffEinsteinRelation() {}

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
