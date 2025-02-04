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
#include "regionProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for external electric field."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting iteration loop\n" << endl;

    // while ( (*std::max_element(voltExtAmpIter.begin(),voltExtAmpIter.end()) || *std::max_element(voltIndIter.begin(),voltIndIter.end())) )
    // {
    Info<< "Time Iteration = " << runTime.timeIndex() << nl << endl;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // Equations for the External Electric Field
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // Reset counters
    regionLoopCounter = 0;
    solverPerformance::debug = 1;
    Info<< "Current Iteration = " << 1 << endl;

    // Loop through different regions to solve the Laplace equations for voltExtAmp
    while (*std::max_element(voltExtAmpIter.begin(), voltExtAmpIter.end()) && regionLoopCounter<100)
    {
        voltExtAmpIter.clear();

        // Solve gas region
        Foam::solverPerformance solvPerfVolt = solve
        (
            fvm::laplacian(voltExtAmp)
        );
        voltExtAmpIter.push_back(solvPerfVolt.nIterations());

        // Solve dielectric regions
        forAll(solidRegions, i)
        {
            #include "setRegionDielectricFields.H"
            Foam::solverPerformance solvPerfVolt = solve 
            (
                fvm::laplacian(voltExtAmp)
            );
            voltExtAmpIter.push_back(solvPerfVolt.nIterations());
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
    Info<< "External Field Region Inner Loops = " << regionLoopCounter << endl;

    // Calculate the electric field EExtAmp
    EExtAmp = -fvc::grad(voltExtAmp);
    forAll(solidRegions, i)
    {
        #include "setRegionDielectricFields.H"
        EExtAmp = -fvc::grad(voltExtAmp);
    }


    // // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // // Equations for the Induced Electric Field
    // // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // // Reset counters
    // regionLoopCounter = 0;
    // solverPerformance::debug = 1;
    // Info<< "Current Iteration = " << 1 << endl;

    // // Loop through different regions to solve the Laplace equations for voltExtAmp
    // while (*std::max_element(voltIndIter.begin(), voltIndIter.end()) && regionLoopCounter<1)
    // {
    //     voltIndIter.clear();

    //     // Solve gas regions
    //     Foam::solverPerformance solvPerfVolt = solve
    //     (
    //         fvm::laplacian(voltInd) + rhoq/e0
    //     );
    //     voltIndIter.push_back(solvPerfVolt.nIterations());

    //     // Solve dielectric regions
    //     forAll(solidRegions, i)
    //     {
    //         #include "setRegionDielectricFields.H"
    //         Foam::solverPerformance solvPerfVolt = solve 
    //         (
    //             fvm::laplacian(voltInd)
    //         );
    //         voltIndIter.push_back(solvPerfVolt.nIterations());
    //     }

    //     // Print performance at custom iteration intervals
    //     regionLoopCounter++;
    //     solverPerformance::debug = 0;
    //     int printPerformance = 500;
    //     if (regionLoopCounter % printPerformance == 0)
    //     {
    //         Info<< "Current Iteration = " << regionLoopCounter << endl;
    //         solverPerformance::debug = 1;
    //     }
    // }
    // Info<< "Induced Field Region Inner Loops = " << regionLoopCounter << endl;

    // // Calculate the electric field EInd
    // EInd = -fvc::grad(voltInd);
    // forAll(solidRegions, i)
    // {
    //     #include "setRegionDielectricFields.H"
    //     EInd = -fvc::grad(voltInd);
    // }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // Calculate for the time varying External Electric Field
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    while( runTime.loop() )
    {

        // Calculate the potential of the external electric field
        // voltExt = voltHandler->calcVoltage();
        voltExt = voltExtAmp*Foam::sin(2*M_PI*(1.0/endTime)*runTime.value());
        forAll(solidRegions, i)
        {
            // voltExtSolid[i] = voltHandlerSolid[i]->calcVoltage();
            #include "setRegionDielectricFields.H"
            voltExt = voltExtAmp*Foam::sin(2*M_PI*(1.0/endTime)*runTime.value());
            // voltExtSolid[i] = Vmax*voltExtAmpSolid[i]*Foam::sin(2*M_PI*freq*runTime.value());
        }

        // Info << "ALL GOOD" << endl;

        // // Calculate the external electric field
        // EExt = -fvc::grad(voltExt);
        // forAll(solidRegions, i)
        // {
        //     #include "setRegionDielectricFields.H"
        //     EExt = -fvc::grad(voltExt);
        // }

        runTime.printExecutionTime(Info);
        runTime.write();
    }


        // // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // // Calculate Local Fields
        // // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // // Calculate local electric potential in all regions
        // volt = voltExtAmp*Foam::sin(2*M_PI*freq*dischargeTime) + voltInd;
        // forAll(solidRegions, i)
        // {
        //     #include "setRegionDielectricFields.H"
        //     volt = voltExtAmp*Foam::sin(2*M_PI*freq*dischargeTime) + voltInd;
        // }
        
        // // Calculate the local electric field
        // E = -fvc::grad(volt);
        // forAll(solidRegions, i)
        // {
        //     #include "setRegionDielectricFields.H"
        //     E = -fvc::grad(volt);
        // }

        // // Calculate Coulomb force
        // Fc = rhoq*E;

        // Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        //     << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        //     << nl << nl << endl;
        runTime.printExecutionTime(Info);
        runTime.write();
    // }

    runTime.writeAndEnd();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
