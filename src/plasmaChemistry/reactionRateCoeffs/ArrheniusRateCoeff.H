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
    Foam::ArrheniusRateCoeff

Description
    This class provides the reads and defines the plasma chemistry.

SourceFiles
    ArrheniusRateCoeff.C

\*---------------------------------------------------------------------------*/

#ifndef ArrheniusRateCoeff_H
#define ArrheniusRateCoeff_H

#include "reactionRateCoeffsBase.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class ArrheniusRateCoeff Declaration
\*---------------------------------------------------------------------------*/

class ArrheniusRateCoeff : public reactionRateCoeffsBase
{
protected:
    scalar A_, B_, C_; 
    word fieldName_;

public:
    TypeName("ArrheniusLaw");

    // Constructor
    ArrheniusRateCoeff(const dictionary& dict, const word& name, plasmaChemistryModel& chemistry)
    {
        read(dict);
    }

    // Member functions
    virtual void read(const dictionary& reactionsSubDict) override;
    virtual void calculate(plasmaChemistryModel& chemistry, const label reactionIndex) const override;
    
    // Virtual destructor
    virtual ~ArrheniusRateCoeff() {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
