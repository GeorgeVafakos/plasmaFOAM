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
    plasmaChemFoam

Description
    Solver streamer plasma chemistry.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createMeshD.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting iteration loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Iteration = " << runTime.timeName() << nl << endl;

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Equations for the Electric Field
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        // // Reset counters
        // counter = 0;
        // solverPerformance::debug = 1;

        // while ((convVoltAext>0 || convVoltDext>0 ) && counter<1)
        // {
            // Air
            Foam::solverPerformance solvPerfVoltA = solve 
            (
                fvm::laplacian(voltA)
            );
            convVoltA = solvPerfVoltA.nIterations();

            // Dielectric
            Foam::solverPerformance solvPerfVoltD = solve
            (
                fvm::laplacian(voltD)
            );
            convVoltD = solvPerfVoltD.nIterations();

        //     counter++;
        //     solverPerformance::debug = 0;
        //     if (counter % 5000 == 0)
        //     {
        //         Info<< "Current Loop = " << counter << endl;
        //         solverPerformance::debug = 1;
        //     }
        // }

        // Info<< "External Field: Region Inner Loops = " << counter << endl;





        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
