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
    Foam::reactionRateCoeffsBase

Description
    Abstract class that provides virstual functions for calculating
    reaction rate coefficients.

\*---------------------------------------------------------------------------*/

#ifndef reactionRateCoeffsBase_H
#define reactionRateCoeffsBase_H

#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Abstract class reactionRateCoeffsBase Declaration
\*---------------------------------------------------------------------------*/

class plasmaChemistryModel;  // forward declaration

class reactionRateCoeffsBase
{
protected:
    //- inline - return the dimension units of volScalarField field
    dimensionedScalar dim(const volScalarField& field) const
    {
        return dimensionedScalar("field", field.dimensions(), scalar(1.0));
    }

public: 
    // Member functions
    //- read the reaction rate parameters from the dictionary
    virtual void read(const dictionary& reactionsSubDict) = 0;

    //- calculate the reaction rate coefficients
    virtual void calculate(plasmaChemistryModel& chemistry) const = 0;

    // Destructor
    virtual ~reactionRateCoeffsBase() {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
