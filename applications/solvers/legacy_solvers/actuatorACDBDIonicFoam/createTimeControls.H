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


// Custom time write controls
bool writeCustomTime =
    runTime.controlDict().lookupOrDefault("writeCustomTime", true);

int printScreenResults =
    runTime.controlDict().lookupOrDefault("printScreenResults", 1);

scalar startTime =
    runTime.controlDict().lookupOrDefault<scalar>("startTime", 0);

scalar endTime =
    runTime.controlDict().lookupOrDefault<scalar>("endTime", 0.0001);

int numWriteTimes =
    runTime.controlDict().lookupOrDefault("numWriteTimes", 15);

scalar P =
    runTime.controlDict().lookupOrDefault<scalar>("P", 0.1);

scalar Q =
    runTime.controlDict().lookupOrDefault<scalar>("Q", 5.0);

int numDischarges =
    runTime.controlDict().lookupOrDefault("numDischarges", 2);

scalar timeBetweenDischarges =
    runTime.controlDict().lookupOrDefault<scalar>("timeBetweenDischarges", 0.5e-5);

// Create normilized count of number of printed times
scalar n[numWriteTimes];
for (int i = 0; i < numWriteTimes; ++i)
    n[i] = 1.0*i/(numWriteTimes-1.0);

// Create stretching function
scalar s[numWriteTimes];
for (int i = 0; i < numWriteTimes; ++i)
    s[i] = P*n[i] + (1.0-P)*(1.0-Foam::tanh(Q*(1.0-n[i]))/Foam::tanh(Q));
    
int multiplier=3;
scalar n1[multiplier*numWriteTimes];
for (int i = 0; i < multiplier*numWriteTimes; ++i)
    n1[i] = 1.0*i/(multiplier*numWriteTimes-1.0);

// Create stretching function
scalar s1[multiplier*numWriteTimes];
for (int i = 0; i < multiplier*numWriteTimes; ++i)
    s1[i] = P*n1[i] + (1.0-P)*(1.0-Foam::tanh(Q*(1.0-n1[i]))/Foam::tanh(Q));

// Counters
int writeTimeCounter = 0;
int dischargeCounter = 1;

// ************************************************************************* //
