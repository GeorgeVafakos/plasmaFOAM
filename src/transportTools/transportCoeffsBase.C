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

#include "transportCoeffsBase.H" 
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(transportCoeffsBase, 0);
defineRunTimeSelectionTable(transportCoeffsBase, mobility);
defineRunTimeSelectionTable(transportCoeffsBase, diffusion);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
transportCoeffsBase::transportCoeffsBase(const fvMesh& mesh, const dictionary& dict, volScalarField& transportCoeff) :
    mesh_(mesh),
    dict_(dict),
    transportCoeff_(transportCoeff),
    E_(mesh.lookupObjectRef<volVectorField>("E")),
    N_(mesh.lookupObjectRef<volScalarField>("NGas")),
    P_(mesh.lookupObjectRef<volScalarField>("pGas")),
    EtoN_("EtoN", 1.0e21*mag(E_/N_)), 
    EtoP_("EtoP", mag(E_/P_)), 
    dimEtoN_("dimEtoN", EtoN_.dimensions(), scalar(1)),
    dimEtoP_("dimEtoP", EtoP_.dimensions(), scalar(1)),
    dimMu_("dimMu", dimensionSet(-1, 0, 2, 0, 0, 1, 0), scalar(1)),
    dimD_("dimD", dimViscosity, scalar(1)),
    dimN_("dimN", N_.dimensions(), scalar(1)),
    dimP_("dimP", P_.dimensions(), scalar(1)),
    isConstant_(false)
{ }


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
autoPtr<transportCoeffsBase> transportCoeffsBase::NewMobility
(
    const fvMesh& mesh,
    const dictionary& dict,
    volScalarField& transportCoeff
)
{
    const word modelType = dict.get<word>("type");

    // Info<< "   Mobility model: " << modelType;

    auto* ctorPtr = mobilityConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup(dict, "transportCoeffsBase (mobility)", modelType, *mobilityConstructorTablePtr_)
            << exit(FatalIOError);
    }

    return autoPtr<transportCoeffsBase>(ctorPtr(mesh, dict, transportCoeff));
}

autoPtr<transportCoeffsBase> transportCoeffsBase::NewDiffusion
(
    const fvMesh& mesh,
    const dictionary& dict,
    volScalarField& transportCoeff
)
{
    const word modelType = dict.get<word>("type");

    // Info<< "    Diffusion model: " << modelType << endl;

    auto* ctorPtr = diffusionConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup(dict, "transportCoeffsBase (diffusion)", modelType, *diffusionConstructorTablePtr_)
            << exit(FatalIOError);
    }

    return autoPtr<transportCoeffsBase>(ctorPtr(mesh, dict, transportCoeff));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

