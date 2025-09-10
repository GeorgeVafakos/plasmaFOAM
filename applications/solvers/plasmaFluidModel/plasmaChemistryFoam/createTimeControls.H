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

Global
    readTimeControls

Description
    Read the control parameters used by addNewDischarge and writeCustomTime

\*---------------------------------------------------------------------------*/

bool adjustTimeStep =
    runTime.controlDict().lookupOrDefault("adjustTimeStep", false);

scalar maxCo =
    runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);

scalar maxDeltaT =
    runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", 0.1);


// Region convergence variables
int regionLoopCounter=0;
std::vector<int> voltEqnIter = {1};

int maxRegionLoops =
    runTime.controlDict().lookupOrDefault("maxRegionLoops", 10000);

int maxRegionLoopsExtAmpEqn =
    runTime.controlDict().lookupOrDefault("maxRegionLoopsExtAmpEqn", 30000);


// Custom time write controls
int printScreenResults =
    runTime.controlDict().lookupOrDefault("printScreenResults", 1);

    
// Custom update of plasma transport and rate coefficients
bool customUpdatePlasmaCoeffs =
    runTime.controlDict().lookupOrDefault("customUpdatePlasmaCoeffs", false);

int plasmaCoeffsUpdateInterval =
    runTime.controlDict().lookupOrDefault("plasmaCoeffsUpdateInterval", 1);


// ************************************************************************* //
