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
    actuatorStreamerACDBDFoam

Description
    Solver for the electrostatic field for a single streamer in an AC-DBD actuator.

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
    #include "createFieldsSolid.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting iteration loop\n" << endl;

    while (runTime.loop() && (convVoltAext>0 || convVoltArho>0 || convVoltDext>0 || convVoltDrho>0 || convVoltIext>0 || convVoltIrho>0) )
    {
        Info<< "Iteration = " << runTime.timeIndex() << nl << endl;

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Equations for the External Electric Field
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        // Reset counters
        conver = 1;
        counter = 0;
        solverPerformance::debug = 1;

        while (conver && counter<30000)
        {
            // Air
            Foam::solverPerformance solvPerfVoltAext = solve 
            (
                fvm::laplacian(voltAext)
            );
            convVoltAext = solvPerfVoltAext.nIterations();

            // Dielectric
            Foam::solverPerformance solvPerfVoltDext = solve
            (
                fvm::laplacian(voltDext)
            );
            convVoltDext = solvPerfVoltDext.nIterations();

            // Insulator
            Foam::solverPerformance solvPerfVoltIext = solve
            (
                fvm::laplacian(voltIext)
            );
            convVoltIext = solvPerfVoltIext.nIterations();

            // Region convergence
            conver = (solvPerfVoltAext.initialResidual()>1.e-6) && (solvPerfVoltDext.initialResidual()>1.e-6) && (solvPerfVoltIext.initialResidual()>1.e-6);
            counter++;
            
            solverPerformance::debug = 0;
            if (counter % 5000 == 0)
            {
                Info<< "Current Loop = " << counter << endl;
                solverPerformance::debug = 1;
            }
        }
        Info<< "External Field: Region Inner Loops = " << counter << endl;



        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Equations for the Induced Electric Field
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        // Reset counters
        conver = 1;
        counter = 0;
        solverPerformance::debug = 1;

        while (conver && counter<30000)
        {
            // Air
            Foam::solverPerformance solvPerfVoltArho = solve 
            (
                fvm::laplacian(voltArho) + rhoq/e0
            );
            convVoltArho = solvPerfVoltArho.nIterations();

            // Dielectric
            Foam::solverPerformance solvPerfVoltDrho = solve
            (
                fvm::laplacian(voltDrho)
            );
            convVoltDrho = solvPerfVoltDrho.nIterations();

            // Insulator
            Foam::solverPerformance solvPerfVoltIrho = solve
            (
                fvm::laplacian(voltIrho)
            );
            convVoltIrho = solvPerfVoltIrho.nIterations();


            // Region convergence
            conver = (solvPerfVoltArho.initialResidual()>1.e-6) && (solvPerfVoltDrho.initialResidual()>1.e-6) && (solvPerfVoltIrho.initialResidual()>1.e-6);
            counter++;
            
            solverPerformance::debug = 0;
            if (counter % 5000 == 0)
            {
                Info<< "Current Loop = " << counter << endl;
                solverPerformance::debug = 1;
            }
        }
        Info<< "Induced Field: Region Inner Loops = " << counter << endl;
        solverPerformance::debug = 1;


        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Calculate Total Fields
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        
        // Calculate total electric potential in all regions
        voltA = voltAext*Foam::sin(2*3.14159*freq*streamerTime) + voltArho;
        voltD = voltDext*Foam::sin(2*3.14159*freq*streamerTime) + voltDrho;
        voltI = voltIext*Foam::sin(2*3.14159*freq*streamerTime) + voltIrho;
        
        // Calculate the electric field
        EAext = -fvc::grad(voltAext);
        EArho = -fvc::grad(voltArho);
        EA    = -fvc::grad(voltA);
        EDext = -fvc::grad(voltDext);
        EDrho = -fvc::grad(voltDrho);
        ED    = -fvc::grad(voltD);
        EIext = -fvc::grad(voltIext);
        EIrho = -fvc::grad(voltIrho);
        EI    = -fvc::grad(voltI);
        Fc    = rhoq*EA;

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
