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
    Foam::constantRateCoeff

Description
    This class provides the reads and defines the plasma chemistry.

SourceFiles
    constantRateCoeff.C

\*---------------------------------------------------------------------------*/

#ifndef constantReactionRateCoeffs_H
#define constantReactionRateCoeffs_H

#include "reactionRateCoeffsBase.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class constantRateCoeff Declaration
\*---------------------------------------------------------------------------*/

class constantRateCoeff : public reactionRateCoeffsBase
{
protected:
    scalar constantValue_;

public:
    TypeName("constantRate");

    // Constructor
    constantRateCoeff(const dictionary& dict, const word& name, plasmaChemistryModel& chemistry);
    // {
    //     read(dict);
    // }

    // Member functions
    // virtual void read(const dictionary& d) override;
    virtual void calculate(plasmaChemistryModel& chemistry, const label reactionIndex) const override;

    
    
    // Virtual destructor
    virtual ~constantRateCoeff() {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
