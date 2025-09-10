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

#ifndef bolsigRateCoeff_H
#define bolsigRateCoeff_H

#include "reactionRateCoeffsBase.H"
#include "interpolation2DTable.H"
#include "tableBounds.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class bolsigRateCoeff Declaration
\*---------------------------------------------------------------------------*/

class bolsigRateCoeff : public reactionRateCoeffsBase
{
protected:
    label reactionIndex_;
    IOdictionary bolsigProperties_;

    PtrList<interpolation2DTable<scalar>> ratesTable_;
    // static interpolation2DTable<scalar>* bolsigTablePtr_;
    // static bool initialized_;

    // Tabulation variables
    const volScalarField& EN_;
    mutable const volScalarField* HePtr_;

public:
    TypeName("bolsig");

    // Constructor
    bolsigRateCoeff
    (
        const dictionary& dict, 
        const word& name, 
        plasmaChemistryModel& chemistry
    );

    // Member functions
    // virtual void read(const dictionary& reactionsSubDict) override;
    virtual void calculate(plasmaChemistryModel& chemistry, const label reactionIndex) const override;
    
    // Virtual destructor
    virtual ~bolsigRateCoeff() {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
