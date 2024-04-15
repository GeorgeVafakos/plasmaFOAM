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
    tubeReactorElectrostatic3RegionsFoam

Description
    Solver for electrostatics in a DBD actuator.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createMeshD.H"
    #include "createMeshI.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting iteration loop\n" << endl;

    while (runTime.loop() && (convVolt>0 || convVoltD>0 || convVoltI>0))
    {
        Info<< "Iteration = " << runTime.timeIndex() << nl << endl;

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Poisson Equations
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        
        // Reset counters
        counter = 0;
        solverPerformance::debug = 1;
        
        while ((convVolt>0 || convVoltD>0 || convVoltI>0) && counter<30000)
        {

            // Inner tube plasma
            Foam::solverPerformance solvPerfVolt = solve 
            (
                fvm::laplacian(volt)// + rhoq/e0
            );
            convVolt = solvPerfVolt.nIterations();

            // Dielectric
            Foam::solverPerformance solvPerfVoltD = solve
            (
                fvm::laplacian(voltD)
            );
            convVoltD = solvPerfVoltD.nIterations();

            // Air
            Foam::solverPerformance solvPerfVoltI = solve
            (
                fvm::laplacian(voltI)
            );
            convVoltI = solvPerfVoltI.nIterations();
            
            counter++;
            solverPerformance::debug = 0;
            if (counter % 100 == 0)
            {
                Info<< "Current Loop = " << counter << endl;
                solverPerformance::debug = 1;
            }
        }

        // Calculate the electric field
        E = -fvc::grad(volt);
        ED = -fvc::grad(voltD);
        EI = -fvc::grad(voltI);

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << nl << endl;
    }

    runTime.writeAndEnd();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
