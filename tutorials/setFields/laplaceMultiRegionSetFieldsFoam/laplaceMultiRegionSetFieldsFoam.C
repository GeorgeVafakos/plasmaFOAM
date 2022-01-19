/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Application
    laplaceMultiRegionSetFieldsFoam

Description
    Solver for Laplace/Poisson equation in two regions.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createMeshA.H"
    #include "createMeshS.H"
    #include "createFields.H"
    #include "createFieldsSolid.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting iteration loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        int conver = 1;
        int counter = 0;
        solverPerformance::debug = 0;

        while (conver )//&& (counter<100))
        {
            Foam::solverPerformance solvPerfVoltA = solve 
            (
                fvm::laplacian(Tair) - S
            );


            Foam::solverPerformance solvPerfVoltD = solve 
            (
                fvm::laplacian(Tsolid) 
            );
            
            conver = (solvPerfVoltA.initialResidual()>1.e-12) || (solvPerfVoltD.initialResidual()>1.e-12);
            counter++;
        }    
        

        // Calculate the heat flux
        qAir = -kAir*fvc::grad(Tair);
        qSolid = -kSolid*fvc::grad(Tsolid);

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << nl << endl;

    }

    Info<< "End\n" << endl;

    runTime.writeAndEnd();

    return 0;
}


// ************************************************************************* //
