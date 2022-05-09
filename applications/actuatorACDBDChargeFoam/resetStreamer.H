/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
     \\/     M anipulation  |
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

if ( runTime.time().value()>( startTime+lengthStreamer*streamerCounter ) && streamerCounter<numStreamers)
{
    rhoq = rhoq_strm;
    voltArho = voltArho_strm;
    voltDrho = voltDrho_strm;
    voltIrho = voltIrho_strm;

    voltA = voltAext*Foam::sin(2*M_PI*(1.0/endTime)*runTime.value()) + voltArho;
    voltD = voltDext*Foam::sin(2*M_PI*(1.0/endTime)*runTime.value()) + voltDrho;
    voltI = voltIext*Foam::sin(2*M_PI*(1.0/endTime)*runTime.value()) + voltIrho;

    EA = -fvc::grad(voltA);
    ED = -fvc::grad(voltD);
    EI = -fvc::grad(voltI);
    Fc = rhoq*EA;

    streamerCounter++;
    writeTimeCounter = 0;
}
