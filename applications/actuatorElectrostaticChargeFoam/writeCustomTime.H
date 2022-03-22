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
    writeCustomTime

Description
    Read the control parameters used by setDeltaT

\*---------------------------------------------------------------------------*/


// char printCustomWriteInterval[]

// startTime

// Info<< "Test strings = " << tempo[2]<< endl;

// // Create normilized count of number of printed times
// for (int i = 0; i < numWriteTimes; ++i)
//     n[i] = i/(numWriteTimes-1.0);

// // Create stretching function
// for (int i = 0; i < numWriteTimes; ++i)
//     s[i] = P*n[i] + (1.0-P)*(1.0-Foam::tanh(Q*(1.0-n[i]))/Foam::tanh(Q));

// system("foamDictionary -entry "writeInterval" -set "0.1e-5" system/controlDict");

// tStr = tA + (P*n+(1-P)*(1-tanh(Q*(1-n))/tanh(Q)))*(tB-tA);


// Info<< "This is the current tiiiiiiiiiime +1 = " << runTime.time().value()+1 << endl;

int i = 1;
if ( runTime.time().value() >= (startTime+s[i]*(endTime-startTime)) )
{
    CmmndStr2 = std::to_string( (s[i]-s[i-1])*(endTime-startTime) );
    system(CmmndStr1+CmmndStr2+CmmndStr3);
    i++;
}

// Info<< str1+str2+str3 << endl;

// CmmndStr2 = std::to_string(0.00012);

// system(CmmndStr1+CmmndStr2+CmmndStr3);


// Info<< "Test s[0] = " << s[0]<< endl;
// Info<< "Test s[1] = " << s[1]<< endl;
// Info<< "Test s[2] = " << s[2]<< endl;
// Info<< "Test s[14] = " << s[14]<< endl;
// Info<< "Test s[14] = " << s[numWriteTimes-1]<< endl;
// Info<< "tanh(1) = " << Foam::tanh(0.1)<< endl;

// system("foamDictionary -entry \"writeInterval\" -set \"0.1e-5\" system/controlDict");
// system("ls");