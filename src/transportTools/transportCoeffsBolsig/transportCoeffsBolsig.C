/*--------------------------# include "transportCoeffsBolsig.H" 
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(transportCoeffsBolsig, 0);
addToRunTimeSelectionTable(transportCoeffsHandler, transportCoeffsBolsig, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
transportCoeffsBolsig::transportCoeffsBolsig(const fvMesh& mesh, const dictionary& dict, volScalarField& transpCoeff):
    transportCoeffsHandler(mesh, dict, transpCoeff)--------------------------------------*\
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

#include "transportCoeffsBolsig.H" 
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(transportCoeffsBolsig, 0);
addToRunTimeSelectionTable(transportCoeffsHandler, transportCoeffsBolsig, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
transportCoeffsBolsig::transportCoeffsBolsig(const fvMesh& mesh, const dictionary& dict, volScalarField& transpCoeff):
    transportCoeffsHandler(mesh, dict, transpCoeff),
    mobilityTable_(),
    diffusivityTable_()
    // mob_e_(transpCoeff),
    // diff_e_(transpCoeff)
{
    // Read the Bolsig properties from the constant directory
    Info<< "Reading bolsigProperties\n" << endl;
    IOdictionary bolsigProperties_
    (
        IOobject
        (
            "bolsigProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    scalarList ENList_
    (
        bolsigProperties_.lookup("EN")
    );

    scalarList HeList_
    (
        bolsigProperties_.lookup("He")
    );

    scalarList eMobList_
    (
        bolsigProperties_.lookup("electronMobility")
    );

    scalarList eDifList_
    (
        bolsigProperties_.lookup("electronDiffusivity")
    );

    // Create lookup tables for mobility and diffusivity
    List<Tuple2<scalar, List<Tuple2<scalar, scalar>>>> mobTable_(HeList_.size());
    List<Tuple2<scalar, List<Tuple2<scalar, scalar>>>> difTable_(HeList_.size());

    forAll(HeList_, i)
    {
        List<Tuple2<scalar, scalar>> mobRow_(ENList_.size());
        List<Tuple2<scalar, scalar>> diffRow_(ENList_.size());

        forAll(ENList_, j)
        {
            label index = i * ENList_.size() + j;
            mobRow_[j]  = Tuple2<scalar, scalar>(ENList_[j], eMobList_[index]);
            diffRow_[j] = Tuple2<scalar, scalar>(ENList_[j], eMobList_[index]);
        }

        mobTable_[i]  = Tuple2<scalar, List<Tuple2<scalar, scalar>>>(HeList_[i], mobRow_);
        difTable_[i] = Tuple2<scalar, List<Tuple2<scalar, scalar>>>(HeList_[i], diffRow_);
    }

    mobilityTable_ = interpolation2DTable<scalar>
    (
        mobTable_,
        bounds::normalBounding(bounds::normalBounding::CLAMP),
        fileName("mobilityTable")
    );

    diffusivityTable_ = interpolation2DTable<scalar>
    (
        difTable_,
        bounds::normalBounding(bounds::normalBounding::CLAMP),
        fileName("diffusivityTable")
    );

    // Info<< "Transport coefficients Bolsig initialized with "
    //     << E_N.size() << " E_N values and "
    //     << He.size() << " He values." << endl;

    Info<< "TEST: Transport coefficients Bolsig initialized." << endl;
}

// defineTypeNameAndDebug(transportCoeffsBolsig, 0);
// addToRunTimeSelectionTable(transportCoeffsHandler, transportCoeffsBolsig, dictionary);

// // 4-arg constructor
// transportCoeffsBolsig::transportCoeffsBolsig
// (
//     const fvMesh& mesh,
//     const dictionary& dict,
//     volScalarField& electronMobilityField,
//     volScalarField& electronDiffusivityField
// )
// :
//     transportCoeffsHandler(mesh, dict, electronMobilityField),
//     mob_e_(electronMobilityField),
//     diff_e_(electronDiffusivityField)
// {
//     // build your mobilityTable_ and diffusivityTable_ here...
// }

// // 3-arg overload (just wraps the 4-arg, giving a dummy for diffusivity)
// transportCoeffsBolsig::transportCoeffsBolsig
// (
//     const fvMesh& mesh,
//     const dictionary& dict,
//     volScalarField& electronMobilityField
// )
// :
//     transportCoeffsBolsig(mesh, dict, electronMobilityField, electronMobilityField)
// {}

// // fill mobility
void transportCoeffsBolsig::calcMobilityCoeffs()
{
    Info<< "TEST: Calculating mobility coefficients." << endl;
    const volScalarField& EN_ = mesh_.lookupObjectRef<volScalarField>("E_N");
    const volScalarField& HeF_ = mesh_.lookupObjectRef<volScalarField>("HeF");

    forAll(transpCoeff_, celli)
    {
        transpCoeff_[celli] = mobilityTable_(HeF_[celli], EN_[celli]);
    }
}

// // fill diffusivity
void transportCoeffsBolsig::calcDiffusionCoeffs()
{
    Info<< "TEST: calcDiffusionCoeffs called." << endl;
//     const volScalarField& EN = lookupVolField<scalar>("E_N");
//     const volScalarField& He = lookupVolField<scalar>("He");

//     forAll(diff_e_, celli)
//     {
//         diff_e_[celli] = diffusivityTable_(EN[celli], He[celli]);
//     }
}




// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// void Foam::transportCoeffsBolsig::calcMobilityCoeffs()
// {
//     const volScalarField& EN = mesh_.lookupObject<volScalarField>("E_N");
//     const volScalarField& He = mesh_.lookupObject<volScalarField>("He");

//     forAll(mob_e_, celli)
//     {
//         mob_e_[celli] = mobilityTable_(EN[celli], He[celli]);
//     }
// }

// void Foam::transportCoeffsBolsig::calcDiffusionCoeffs()
// {
//     const volScalarField& EN = mesh_.lookupObject<volScalarField>("E_N");
//     const volScalarField& He = mesh_.lookupObject<volScalarField>("He");

//     forAll(diff_e_, celli)
//     {
//         diff_e_[celli] = diffusivityTable_(EN[celli], He[celli]);
//     }
// }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



