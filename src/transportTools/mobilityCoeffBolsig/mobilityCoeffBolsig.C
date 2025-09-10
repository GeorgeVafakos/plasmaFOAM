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

#include "mobilityCoeffBolsig.H" 
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(mobilityCoeffBolsig, 0);
addToRunTimeSelectionTable(transportCoeffsBase, mobilityCoeffBolsig, mobility);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
mobilityCoeffBolsig::mobilityCoeffBolsig(const fvMesh& mesh, const dictionary& dict, volScalarField& transportCoeff):
    transportCoeffsBase(mesh, dict, transportCoeff),
    mobilityTable_(),
    EN_(mesh_.lookupObjectRef<volScalarField>("EN")),
    He_(mesh_.lookupObjectRef<volScalarField>("He"))
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

    // Create lookup tables for mobility and diffusivity
    List<Tuple2<scalar, List<Tuple2<scalar, scalar>>>> mobTable_(HeList_.size());

    forAll(HeList_, i)
    {
        List<Tuple2<scalar, scalar>> mobRow_(ENList_.size());

        forAll(ENList_, j)
        {
            label index = i * ENList_.size() + j;
            mobRow_[j]  = Tuple2<scalar, scalar>(ENList_[j], eMobList_[index]);
        }

        mobTable_[i]  = Tuple2<scalar, List<Tuple2<scalar, scalar>>>(HeList_[i], mobRow_);
    }

    mobilityTable_ = interpolation2DTable<scalar>
    (
        mobTable_,
        bounds::normalBounding(bounds::normalBounding::CLAMP),
        fileName("mobilityTable")
    );
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void mobilityCoeffBolsig::calcMobilityCoeffs()
{
    forAll(transportCoeff_, celli)
    {
        transportCoeff_[celli] = mobilityTable_(He_[celli], EN_[celli]);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



