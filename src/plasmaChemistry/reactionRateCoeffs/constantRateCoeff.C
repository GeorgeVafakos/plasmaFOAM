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

\*---------------------------------------------------------------------------*/

#include "constantRateCoeff.H"
#include "plasmaChemistryModel.H"
#include "volFields.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

void constantRateCoeff::calculate(plasmaChemistryModel& chemistry) const
{
    PtrList<volScalarField>& k = chemistry.k();
    const IOdictionary& reactionsDict = chemistry.reactionsDict();
    const dictionary& reactions(reactionsDict.subDict("reactions"));

    label j = 0;
    for (const entry& e : reactions)
    {
        const dictionary& reacDict = e.dict();
        // const dictionary& reacDict = iter().dict();

        word type = word(reacDict.lookup("type"));

        if (type != "constantRate")
        {
            ++j;
            continue;  // Skip if not this type
        }

        // Read the constant value or use default
        scalar coeffValue = reacDict.lookupOrDefault<scalar>("value", 1e-10);

        // Assign value to internalField of k[j]
        k[j] = dim(k[j]) * coeffValue;
        // k[j].correctBoundaryConditions();

        Info<< "Set constant k[" << j << "] = " << coeffValue
            << " for reaction " << e.keyword() << " " << k[j].dimensions() << nl;

        ++j;
    }
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

