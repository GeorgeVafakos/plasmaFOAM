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

Local
    addNewDischarge

Description
    Superimpose additional discharges over existing electrostatic conditions

\*---------------------------------------------------------------------------*/

if (runTime.time().value() >= (startTime+timeBetweenDischarges*dischargeCounter) && dischargeCounter<numDischarges)
{
    // Gas region
    rhoq = rhoq+rhoqDischarge;
    voltInd = voltInd+voltIndDischarge;
    voltExt = voltExtAmp*Foam::sin(2*M_PI*(1.0/endTime)*runTime.value());
    volt = voltExt + voltInd;
    EExt = -fvc::grad(voltExt);
    EInd = -fvc::grad(voltInd);
    E = -fvc::grad(volt);
    Fc = rhoq*E;

    // Dielectric regions
    forAll(solidRegions, i)
    {
        #include "setRegionDielectricFields.H"
        voltInd = voltInd+voltIndDischarge;
        voltExt = voltExtAmp*Foam::sin(2*M_PI*(1.0/endTime)*runTime.value());
        volt = voltExt + voltInd;
        EExt = -fvc::grad(voltExt);
        EInd = -fvc::grad(voltInd);
        E = -fvc::grad(volt);
    }

    dischargeCounter++;
    writeTimeCounter = 0;
    runTime.writeNow();
}
