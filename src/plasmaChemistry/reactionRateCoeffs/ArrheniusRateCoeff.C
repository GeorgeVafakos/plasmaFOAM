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

#include "ArrheniusRateCoeff.H"
#include "plasmaChemistryModel.H"
#include "volFields.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

void ArrheniusRateCoeff::calculate(plasmaChemistryModel& chemistry) const
{
    PtrList<volScalarField>& k = chemistry.k();
    const IOdictionary& reactionsDict = chemistry.reactionsDict();
    const dictionary& reactionsList(reactionsDict.subDict("reactions"));

    label j = 0;
    for (const entry& e : reactionsList)
    {
        const dictionary& reaction = e.dict();
        
        word type = word(reaction.lookup("type"));
        if (type != "ArrheniusLaw")
        {
            ++j;
            continue;  // Skip if not this type
        }

        // Get the field variable name
        word Tword = word(reaction.lookup("field"));

        // Get temperature from the mesh
        const volScalarField& T = chemistry.mesh().lookupObjectRef<volScalarField>(Tword);

        // Read the constant value or use default
        scalar A = reaction.lookupOrDefault<scalar>("A", 1e-10);
        scalar B = reaction.lookupOrDefault<scalar>("B", 1e-10);
        scalar C = reaction.lookupOrDefault<scalar>("C", 1e-10);

        // Assign value to internalField of k[j]
        k[j] = dim(k[j]) * A*pow(T/dim(T), B)*exp(-C / (T/dim(T)));

        // Info<< "DIM T = " << dim(T) << endl;
        // k[j].correctBoundaryConditions();

        // Info<< "Set constant k[" << j << "] = " << coeffValue
        //     << " for reaction " << e.keyword() << " " << k[j].dimensions() << nl;

        ++j;
    }
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

