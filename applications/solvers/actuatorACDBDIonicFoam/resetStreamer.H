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
    resetStreamer

Description
    Initialize new streamer electrostatic conditions

\*---------------------------------------------------------------------------*/

if ( runTime.time().value()>=( startTime+lengthStreamer*streamerCounter ) && streamerCounter<numStreamers)
{
    rhoq = rhoq+rhoq_strm;
    voltArho = voltArho+voltArho_strm;
    voltDrho = voltDrho+voltDrho_strm;
    voltErho = voltErho+voltErho_strm;

    voltAext = voltAextMag*Foam::sin(2*M_PI*(1.0/endTime)*runTime.value());
    voltDext = voltDextMag*Foam::sin(2*M_PI*(1.0/endTime)*runTime.value());
    voltEext = voltEextMag*Foam::sin(2*M_PI*(1.0/endTime)*runTime.value());
        
    voltA = voltAext + voltArho;
    voltD = voltDext + voltDrho;
    voltE = voltEext + voltErho;

    EAext = -fvc::grad(voltAext);
    EArho = -fvc::grad(voltArho);
    EA    = -fvc::grad(voltA);
    EDext = -fvc::grad(voltDext);
    EDrho = -fvc::grad(voltDrho);
    ED    = -fvc::grad(voltD);
    EEext = -fvc::grad(voltEext);
    EErho = -fvc::grad(voltErho);
    EE    = -fvc::grad(voltE);
    Fc = rhoq*EA;

    streamerCounter++;
    writeTimeCounter = 0;
    runTime.writeNow();
}
