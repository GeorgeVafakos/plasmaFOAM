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
    inertSpecie_(reactionsDict_.lookup("inertSpecie")),
    k_(reactions_.size())
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

    // Create stoichiometric matrices
    buildStoichiometricMatrices();

    // Define the reaction rate coefficients
    forAll(reactions_, j)
    {
        k_.set
        (
            j,
            new volScalarField
            (
                IOobject
                (
                    "k" + name(j),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", dimensionSet(0, -3*(sumReactants(j)-1), -1, 0, 0, 0, 0), scalar(1.0e-10))
            )
        );
    }



    label i = 0;
    forAllConstIter(dictionary, reactionsSubDict, iter)
    {
        const dictionary& d = iter().dict();
        word type = d.lookup("type");

        autoPtr<reactionRateCoeffsBase> rateModel;

        if (type == "constantRate")
        {
            rateModel.reset(new constantRateCoeff());
        }
        else if (type == "ArrheniusLaw")
        {
            rateModel.reset(new ArrheniusRateCoeff());
        }
        else
        {
            FatalErrorIn("plasmaChemistryModel") << "Unknown reaction type: " << type << abort(FatalError);
        }

        rateModel().read(d);
        rateCalculators_.set(i++, std::move(rateModel));
}

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

scalar plasmaChemistryModel::sumReactants(const label j) const
{
    scalar sum = 0.0;
    forAll(aL_, i)
    {
        sum += aL_[i][j];
    }
    return sum;
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


volScalarField plasmaChemistryModel::R(const label i, const PtrList<volScalarField>& n) const
{
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

    forAll(reactions_, j)
    {
        // Initialize rate as the rate coefficient
        // scalarField rate(k_[j]); // this also works, must check why
        scalarField rate(k_[j].internalField());

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

    // TO DO: See if i can remove .correctBoundaryConditions() and .internalFieldRef().field()

    return Ri;
}







// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
