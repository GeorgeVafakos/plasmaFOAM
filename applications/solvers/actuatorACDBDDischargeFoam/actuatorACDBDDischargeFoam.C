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
    actuatorACDBDDischargeFoam

Description
    Solver for the electrostatic field for a single discharge event for 
    an AC-DBD actuator.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createMeshD.H"
    #include "createMeshE.H"
    #include "createFields.H"
    #include "createFieldsSolid.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting iteration loop\n" << endl;

    while (runTime.loop() && (convVoltAextMag>0 || convVoltArho>0 || convVoltDextMag>0 || convVoltDrho>0 || convVoltEextMag>0 || convVoltErho>0) )
    {
        Info<< "Iteration = " << runTime.timeIndex() << nl << endl;

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Equations for the External Electric Field
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        // Reset counters
        counter = 0;
        solverPerformance::debug = 1;

        while ((convVoltAextMag>0 || convVoltDextMag>0 || convVoltEextMag>0) && counter<30000)
        {
            // Air
            Foam::solverPerformance solvPerfVoltAextMag = solve 
            (
                fvm::laplacian(voltAextMag)
            );
            convVoltAextMag = solvPerfVoltAextMag.nIterations();

            // Dielectric
            Foam::solverPerformance solvPerfVoltDextMag = solve
            (
                fvm::laplacian(voltDextMag)
            );
            convVoltDextMag = solvPerfVoltDextMag.nIterations();

            // Encapsulator
            Foam::solverPerformance solvPerfVoltEextMag = solve
            (
                fvm::laplacian(voltEextMag)
            );
            convVoltEextMag = solvPerfVoltEextMag.nIterations();

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
        counter = 0;
        solverPerformance::debug = 1;

        while ((convVoltArho>0 || convVoltDrho>0 || convVoltErho>0) && counter<30000)
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

            // Encapsulator
            Foam::solverPerformance solvPerfVoltErho = solve
            (
                fvm::laplacian(voltErho)
            );
            convVoltErho = solvPerfVoltErho.nIterations();

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
        voltA = voltAextMag*Foam::sin(2*3.14159*freq*dischargeTime) + voltArho;
        voltD = voltDextMag*Foam::sin(2*3.14159*freq*dischargeTime) + voltDrho;
        voltE = voltEextMag*Foam::sin(2*3.14159*freq*dischargeTime) + voltErho;
        
        // Calculate the electric field
        EAextMag = -fvc::grad(voltAextMag);
        EArho    = -fvc::grad(voltArho);
        EA       = -fvc::grad(voltA);
        EDextMag = -fvc::grad(voltDextMag);
        EDrho    = -fvc::grad(voltDrho);
        ED       = -fvc::grad(voltD);
        EEextMag = -fvc::grad(voltEextMag);
        EErho    = -fvc::grad(voltErho);
        EE       = -fvc::grad(voltE);
        Fc       = rhoq*EA;

        // runTime.write();

        // Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        //     << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        //     << nl << nl << endl;
        runTime.printExecutionTime(Info);
        runTime.write();
    }

    runTime.writeAndEnd();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
