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

#include "plasmaChemistryModel.H"
#include "Time.H"
#include "volFields.H"
#include "OStringStream.H"
#include "HashPtrTable.H"
#include "speciesTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
    const word Reaction<dummyThermo>::typeName("Reaction<dummyThermo>");

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
plasmaChemistryModel::plasmaChemistryModel(const fvMesh& mesh)
:
    mesh_(mesh),
    reactionsDict_
    (
        IOobject
        (
            "reactionsDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    species_(reactionsDict_.lookup("species")),
    reactions_(reactionsDict_.subDict("reactions").size()),
    inertSpecie_(reactionsDict_.lookup("inertSpecie"))
{
    // Check if inert specie exists
    if (!species_.found(inertSpecie_))
    {
        FatalIOErrorIn("plasmaChemistryModel", reactionsDict_)
            << "Inert specie '" << inertSpecie_ << "' not found in available species: " 
            << species_ << exit(FatalIOError);
    }

    // Read the reactions from the dictionary
    Info << "Reading reactionsDict" << endl;
    const dictionary& reactionsSubDict = reactionsDict_.subDict("reactions");
    speciesHash_ = hashedWordList(species_, true);

    label i = 0;
    forAllConstIter(dictionary, reactionsSubDict, iter)
    {
        HashPtrTable<dummyThermo> dummyTable;

        forAll(species_, j)
        {
            dummyTable.insert(species_[j], autoPtr<dummyThermo>(new dummyThermo(species_[j])));
        }

        reactions_.set
        (
            i,
            new Reaction<dummyThermo>
            (
                speciesHash_,
                dummyTable,
                iter().dict(),
                false,
                false
            )
        );
        ++i;
    }

    Info<< "---------- HELlOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO 2----------" << endl;

    // Create stoichiometric matrices
    buildStoichiometricMatrices();

}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void plasmaChemistryModel::buildStoichiometricMatrices()
{
    const label nSpec = species_.size();
    const label nReac = reactions_.size();

    aL_.setSize(nSpec);
    aR_.setSize(nSpec);
    stoichMatrix_.setSize(nSpec);
    reactantMap_.setSize(nReac);

    forAll(aL_, i)
    {
        aL_[i].setSize(nReac, 0.0);
        aR_[i].setSize(nReac, 0.0);
        stoichMatrix_[i].setSize(nReac, 0.0);
    }

    forAll(reactions_, j)
    {
        const Reaction<dummyThermo>& reaction = reactions_[j];

        List<Tuple2<label, scalar>> nonZeroReactants;

        // Construct aL
        for (const Reaction<dummyThermo>::specieCoeffs& r : reaction.lhs())
        {
            aL_[r.index][j] = r.stoichCoeff;

            if (r.stoichCoeff > SMALL)
            {
                nonZeroReactants.append(Tuple2<label, scalar>(r.index, r.stoichCoeff));
            }
        }
        reactantMap_[j] = nonZeroReactants;

        // Construct aR
        for (const Reaction<dummyThermo>::specieCoeffs& p : reaction.rhs())
        {
            aR_[p.index][j] = p.stoichCoeff;
        }

    }

    // Compute net stoichiometry S
    forAll(stoichMatrix_, i)
    {
        stoichMatrix_[i] = aR_[i] - aL_[i];

    }

}

void plasmaChemistryModel::printStoichiometricMatrices() const
{
    Info << "aL (reactants):" << nl;
    forAll(aL_, i)
    {
        forAll(aL_[i], j)
        {
            Info << int(aL_[i][j]) << " ";
        }
        Info << "// " << species_[i] << nl;
    }

    Info << "aR (products):" << nl;
    forAll(aR_, i)
    {
        forAll(aR_[i], j)
        {
            Info << int(aR_[i][j]) << " ";
        }
        Info << "// " << species_[i] << nl;
    }

    Info << "S (net = aR - aL):" << nl;
    forAll(stoichMatrix_, i)
    {
        forAll(stoichMatrix_[i], j)
        {
            Info << int(stoichMatrix_[i][j]) << " ";
        }
        Info << "// " << species_[i] << nl;
    }
}



// PtrList<volScalarField> plasmaChemistryModel::getSpeciesFields() const
// {
//     PtrList<volScalarField> n(species_.size());
// 
//     forAll(species_, i)
//     {
//         const word& name = species_[i];
// 
//         if (mesh_.foundObject<volScalarField>(name))
//         {
//             const volScalarField& field = mesh_.lookupObject<volScalarField>(name);
//             n.set(i, const_cast<volScalarField*>(&field));
//         }
//         else
//         {
//             FatalErrorInFunction
//                 << "Species field '" << name << "' not found in mesh registry."
//                 << exit(FatalError);
//         }
//     }
// 
//     return n;
// }


volScalarField plasmaChemistryModel::R(const label i, const PtrList<volScalarField>& n) const
{
    // const fvMesh& mesh = mesh_;
    const label nCells = mesh_.nCells();
    const label nReactions = reactions_.size();
    // const label nSpecies = species_.size();

    // Create result field (initialized to 0)
    volScalarField Ri
    (
        IOobject
        (
            "R_" + species_[i],
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimensionSet(0,-3,-1,0,0,0,0), 0.0)
    );

    const scalar kf = 1.0e-10;  // Example constant rate coefficient

    for (label j = 0; j < nReactions; ++j)
    {
        // Initialize rate as the rate coefficient
        scalarField rate(nCells, kf);

        // Get the non-zero reactants for reaction j
        const List<Tuple2<label, scalar>>& nonZeroReactants = reactantMap_[j];

        // Calculate rate as the product of rate coefficient times the reactants raised to their stoichiometric coefficients

        // for (label l = 0; l < nSpecies; ++l)
        // {
        //     const scalar nuL = aL_[l][j];
        //     if (nuL > SMALL)
        //     {
        //         const scalarField& nl = n[l].internalField();
        //         rate *= pow(nl, nuL);
        //     }
        // }

        forAll(nonZeroReactants, k)
        {
            const label l = nonZeroReactants[k].first();
            const scalar nuL = nonZeroReactants[k].second();
            rate *= pow(n[l].internalField(), nuL);
        }

        // Multiply the rate with the stoichiometric coefficient for species i
        scalar Sij = stoichMatrix_[i][j];
        if (mag(Sij) > SMALL)
        {
            Ri.internalFieldRef().field() += Sij * rate;
        }
    }

    Ri.correctBoundaryConditions();
    return Ri;
}













// PtrList<volScalarField> plasmaChemistryModel::R
// (
//     const PtrList<volScalarField>& ni
// ) const
// {
//     const label nSpecies = species_.size();
//     const label nReactions = reactions_.size();
//     const label nCells = mesh_.nCells();

//     // Create result list for species source terms
//     PtrList<volScalarField> Ri(nSpecies);

//     // Placeholder constant rate coefficient
//     const scalar kf = 1.0e-10;

//     // Precompute forward rate term (product of reactants raised to stoich coeffs)
//     List<scalarField> reactionRates(nReactions);

//     for (label r = 0; r < nReactions; ++r)
//     {
//         scalarField Rr(nCells, kf);

//         for (label j = 0; j < nSpecies; ++j)
//         {
//             const scalar nuL = aL_[j][r];
//             if (nuL > SMALL)
//             {
//                 const scalarField& nj = ni[j].internalField();
//                 Rr *= Foam::pow(nj, nuL);
//             }
//         }

//         reactionRates[r] = Rr;
//     }

//     // Build total source term per species using stoichiometry
//     for (label i = 0; i < nSpecies; ++i)
//     {
//         scalarField RiField(nCells, 0.0);

//         for (label r = 0; r < nReactions; ++r)
//         {
//             const scalar S_ir = stoichMatrix_[i][r];
//             if (mag(S_ir) > SMALL)
//             {
//                 RiField += S_ir * reactionRates[r];
//             }
//         }

//         // Construct volScalarField with appropriate dimensions and patch types
//         Ri.set(i, new volScalarField
//         (
//             IOobject
//             (
//                 "R_" + species_[i],
//                 mesh_.time().timeName(),
//                 mesh_,
//                 IOobject::NO_READ,
//                 IOobject::NO_WRITE
//             ),
//             mesh_,
//             dimensionedScalar("zero", dimMass/dimTime, 0.0),  // Adjust if needed
//             ni[i].boundaryField().types()
//         ));

//         Ri[i].internalFieldRef().field() = RiField;
//         Ri[i].correctBoundaryConditions();
//     }

//     return Ri;
// }










// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
