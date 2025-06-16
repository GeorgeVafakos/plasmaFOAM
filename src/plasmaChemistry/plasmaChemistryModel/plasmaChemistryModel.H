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
    Foam::plasmaChemistryModel

Description
    This class provides the reads and defines the plasma chemistry.

SourceFiles
    plasmaChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef plasmaChemistryModel_H
#define plasmaChemistryModel_H

#include "fvMesh.H"
#include "Reaction.H"
#include "dummyThermo.H"
#include "volFieldsFwd.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class plasmaChemistryModel Declaration
\*---------------------------------------------------------------------------*/

class plasmaChemistryModel
{
protected:
    const fvMesh& mesh_;
    IOdictionary reactionsDict_;
    wordList species_;
    hashedWordList speciesHash_;
    PtrList<Reaction<dummyThermo>> reactions_;
    word inertSpecie_;
    List<scalarField> aL_, aR_, stoichMatrix_;
    List<List<Tuple2<label, scalar>>> reactantMap_;

public:

    // Constructor
    plasmaChemistryModel(const fvMesh& mesh);

    // Member operators

    //- inline - get the species 
    const wordList& species() const
    {
        return species_;
    }

    //- inline - get the inert species index
    int inertIndex() const
    {
        return species_.find(inertSpecie_);
    }

    //- inline - get the reacitons
    void printReactions() const
    {
        forAll(reactions_, i)
        {
            Info << "Reaction " << i << ": " << reactions_[i];
        }
    }

    //- print stoichiometric matrices
    void printStoichiometricMatrices() const;

    //- create the stoichiometric matrix
    void buildStoichiometricMatrices();

    // //- Reconstruct list of species fields from the mesh
    // PtrList<volScalarField> getSpeciesFields() const;

    // Return total reaction source term for species with index iSpecie
    volScalarField R(const label iSpecie, const PtrList<volScalarField>& n) const;
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
