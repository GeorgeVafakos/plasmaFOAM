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
    Foam::transportCoeffsHeIonLXCat

Description
    This is a class that provides the calculation of the transport
    coefficients (mobility and diffusion) for the ions.

SourceFiles
    transportCoeffsHeIonLXCat.C

\*---------------------------------------------------------------------------*/

#ifndef transportCoeffsHeIonLXCat_H
#define transportCoeffsHeIonLXCat_H

#include "transportCoeffsBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class transportCoeffsHeIonLXCat Declaration
\*---------------------------------------------------------------------------*/

class transportCoeffsHeIonLXCat : public transportCoeffsBase
{
protected:  
    // Name of the E/N field
    word fieldName_;

    // Reduced electric field (Td)
    volScalarField& EN_;

public:
    TypeName("HeIonLXCat");

    // Constructor
    transportCoeffsHeIonLXCat(const fvMesh& mesh, const dictionary& dict, volScalarField& transportCoeff);

    // Destructor
    ~transportCoeffsHeIonLXCat() {}

    // member functions
    //- calculate mobility coefficients
    void calcMobilityCoeffs();

    //- calculate diffusion coefficients
    void calcDiffusionCoeffs();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
