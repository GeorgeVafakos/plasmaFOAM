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
    Foam::bolsigRateCoeff

Description
    This class provides the reads and defines the plasma chemistry.

SourceFiles
    bolsigRateCoeff.C

\*---------------------------------------------------------------------------*/

#ifndef bolsig3BodyRateCoeff_H
#define bolsig3BodyRateCoeff_H

#include "reactionRateCoeffsBase.H"
#include "bolsigRateCoeff.H"

// #include "interpolation2DTable.H"
// #include "tableBounds.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class bolsigRateCoeff Declaration
\*---------------------------------------------------------------------------*/

class bolsig3BodyRateCoeff : public bolsigRateCoeff
{
protected:
    const volScalarField& NGas_;

public:
    TypeName("bolsig3BodyReaction");

    // Constructor
    bolsig3BodyRateCoeff
    (
        const dictionary& dict, 
        const word& name, 
        plasmaChemistryModel& chemistry
    );

    virtual ~bolsig3BodyRateCoeff() {}

    virtual void calculate(plasmaChemistryModel& chemistry, const label reactionIndex) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
