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

Application
    threeSpeciesSDBDFoam

Description
    Solver for the plasma discharge using the three species model.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "regionProperties.H"
#include "voltageHandler.H"
//#include "transportCoeffsHandler.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for the local electric field."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "createVoltExtObjs.H"

    // Calculate amplitude of external field
    #include "EExtAmpEqn.H"

    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;
    while( runTime.loop() )
    {
        // Calculate the induced electric field
        #include "EIndEqn.H"
        
        // Calculate the local electric field
        #include "calcE.H"

        // Control time step according to Co num
        #include "CourantNo.H"
        #include "setDeltaT.H" 

        ++runTime;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Equations for species continuity        
        #include "NEqn.H"

        runTime.printExecutionTime(Info);
        runTime.write();
    }

    runTime.writeAndEnd();

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
