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
    electricFieldFoam

Description
    Solver for the local electric field.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "regionProperties.H"
#include "voltageHandler.H"


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
    #include "createFields.H"
    #include "createVoltExtObjs.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // Equations for the External Electric Field Amplitude
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting external electric field loop\n" << endl;

    // Reset counters
    regionLoopCounter = 0;
    solverPerformance::debug = 1;
    Info<< "Current Iteration = " << 1 << endl;

    // Loop through different regions to solve the Laplace equations for voltExtAmp
    while (*std::max_element(voltEqnIter.begin(), voltEqnIter.end()) && regionLoopCounter<10)
    {
        voltEqnIter.clear();

        // Solve gas region
        Foam::solverPerformance solvPerfVoltEqn = solve
        (
            fvm::laplacian(voltExtAmp)
        );
        voltEqnIter.push_back(solvPerfVoltEqn.nIterations());

        // Solve solid regions
        forAll(solidRegions, i)
        {
            #include "setRegionSolidFields.H"
            Foam::solverPerformance solvPerfVoltEqn = solve 
            (
                fvm::laplacian(voltExtAmp)
            );
            voltEqnIter.push_back(solvPerfVoltEqn.nIterations());
        }

        // Print performance at custom iteration intervals
        regionLoopCounter++;
        solverPerformance::debug = 0;
        int printPerformance = 500;
        if (regionLoopCounter % printPerformance == 0)
        {
            Info<< "Current Iteration = " << regionLoopCounter << endl;
            solverPerformance::debug = 1;
        }
    }
    Info<< "External Field Amplitude: Region Inner Loops = " << regionLoopCounter << endl;

    // Calculate the external electric field amplitude
    EExtAmp = -fvc::grad(voltExtAmp);
    forAll(solidRegions, i)
    {
        #include "setRegionSolidFields.H"
        EExtAmp = -fvc::grad(voltExtAmp);
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    while( runTime.loop() )
    {
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Calculate for the time varying External Electric Field
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Calculate the potential of the external electric field
        voltExt = voltHandler->calcVoltage();
        forAll(solidRegions, i)
        {
            #include "setRegionSolidFields.H"
            voltExtSolid[i] = voltHandlerSolid[i].calcVoltage();
        }

        // Calculate the external electric field
        EExt = -fvc::grad(voltExt);
        forAll(solidRegions, i)
        {
            #include "setRegionSolidFields.H"
            EExt = -fvc::grad(voltExt);
        }

<<<<<<< HEAD:applications/solvers/externalFieldFoam/externalFieldFoam.C
        // Calculate local fields
=======
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Equations for the time varying Induced Electric Field 
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Reset counters
        regionLoopCounter = 0;
        solverPerformance::debug = 1;
        Info<< "Current Iteration = " << 1 << endl;

        // Loop through different regions to solve the Laplace equations for voltInd
        while (*std::max_element(voltEqnIter.begin(), voltEqnIter.end()) && regionLoopCounter<10)
        {
            voltEqnIter.clear();

            // Solve gas region
            Foam::solverPerformance solvPerfVolt = solve
            (
                fvm::laplacian(voltInd) + rhoq/e0
            );
            voltEqnIter.push_back(solvPerfVolt.nIterations());

            // Solve solid regions
            forAll(solidRegions, i)
            {
                #include "setRegionSolidFields.H"
                Foam::solverPerformance solvPerfVoltEqn = solve 
                (
                    fvm::laplacian(voltInd)
                );
                voltEqnIter.push_back(solvPerfVoltEqn.nIterations());
            }

            // Print performance at custom iteration intervals
            regionLoopCounter++;
            solverPerformance::debug = 0;
            int printPerformance = 500;
            if (regionLoopCounter % printPerformance == 0)
            {
                Info<< "Current Iteration = " << regionLoopCounter << endl;
                solverPerformance::debug = 1;
            }
        }
        Info<< "Induced Field: Region Inner Loops = " << regionLoopCounter << endl;

        // Calculate the external electric field amplitude
        EInd = -fvc::grad(voltInd);
        forAll(solidRegions, i)
        {
            #include "setRegionSolidFields.H"
            EInd = -fvc::grad(voltInd);
        }

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Total fields
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        volt = voltExt + voltInd;
        E = EExt + EInd;
        forAll(solidRegions, i)
        {
            #include "setRegionSolidFields.H"
            volt = voltExt + voltInd;
            E = EExt + EInd;
        }

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Equations for species continuity
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
>>>>>>> driftDiffusionEqn:applications/solvers/electricFieldFoam/electricFieldFoam.C


        runTime.printExecutionTime(Info);
        runTime.write();
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //